# To replicate the simu results in sec 4 of the paper, with known sigma2.

rm(list=ls())
library(intervals)
library(nleqslv)
library(foreach)
library(leaps)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)
library(bigstatsr)

#==============PART 1: helper functions===============
# perform post-selection correction for confidence inference
postICci=function(X_dt, y, selected, alls, new_xpoint, criteria='aic', alpha=0.05, sigma2){
  #input:
  #- X_dt: design matrix
  #- y: response
  #- selected: the AIC selected model (vector: boolean element indicating whether the predictor is included or not)
  #- alls: all other models to be compard with AIC-selected model (matrix: each row represent a model)
  #- new_xpoint: new x point
  #- criteria: model selection criteria, must be one of 'aic', 'bic', 'aicc'
  #- alpha: significance level
  #- sigma2: ground truth noise level in the model
  
  #output:
  #- the post-selection confidence interval for predicted mean at new_xpoint
  
  n=length(y)
  if(criteria=='aic'){
    cons=2
  }else if(criteria=='bic'){
    cons=log(n)
  }else if(criteria=='aicc'){
    #pass
  }
  else{
    print("must be either aic or bic!")
    return(NA)
  }
  #record the excluded intervals for those TruncNormal with 2 intervals
  exclude_lb=rep(NA,nrow(alls))
  exclude_ub=exclude_lb
  #record the included intervals for those TruncNormal with 1 interval
  lb_max=-Inf
  ub_min=Inf
  
  eta=cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%c(1,new_xpoint[selected])
  cee=eta/c(t(eta)%*%eta)
  z=y-cee%*%t(eta)%*%y
  P_select=diag(n)-cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
  
  for (s in 1:nrow(alls)) {
    #apply theorem 2, skip comparing with the selected model itself or its superset
    if(all(alls[s,]-selected>=0)){
      next
    }
    
    #exp weight involving cons
    s_hat=sum(selected)#size of selected model
    s_tp=sum(alls[s,])#size of model alls[s,]
    wt=exp(cons*(s_hat-s_tp)/n)
    if(criteria=='aicc'){
      wt=exp((2*(s_hat-s_tp)+2*s_hat*(s_hat+1)/(n-s_hat-1)-2*s_tp*(s_tp+1)/(n-s_tp-1))/n)
    }
    #P matrix for model s
    P_s=diag(n)-cbind(1,X_dt[,alls[s,]])%*%solve(t(cbind(1,X_dt[,alls[s,]]))%*%cbind(1,X_dt[,alls[s,]]))%*%t(cbind(1,X_dt[,alls[s,]]))
    #find the coef in quadratic form: A*x^2+B*x+C>0
    Pdiff=P_s-P_select*wt
    A=c(t(cee)%*%Pdiff%*%cee)
    B1=c(t(z)%*%Pdiff%*%cee*2)
    C=c(t(z)%*%Pdiff%*%z)
    #solve the quadratic inequality
    Delta=B1^2 - 4 * A * C
    if(Delta>0){
      r1=(-B1 - sqrt(Delta)) / (2 * A)
      r2=(-B1 + sqrt(Delta)) / (2 * A)
      if(A>0){
        #then there must be 2 intervals as domains, and they have one of {-Inf, Inf} as the tip end
        exclude_lb[s]=min(r1,r2)
        exclude_ub[s]=max(r1,r2)
      }else{
        lb_max=max(min(r1,r2),lb_max)
        ub_min=min(max(r1,r2),ub_min)
        
      }
    }else{
      if(A<0){
        lb_max=0
        ub_min=0
        return(NA)
      }#otherwise, the interval is the whole R line.
    }
    
  }
  #remove NAs so that what's left always come in pairs
  #the domain is the intersection of {(lb[2i-1], ub[2i-1])union(lb[2i], ub[2i])}
  exclude_lb=exclude_lb[!is.na(exclude_lb)]
  exclude_ub=exclude_ub[!is.na(exclude_ub)]
  #excluded intervals
  exclude_intervals=interval_union(Intervals(cbind(exclude_lb,exclude_ub)))
  #also note intervals outside (lb_max, ub_min) are also excluded
  if(lb_max>-Inf){
    exclude_intervals=interval_union(exclude_intervals,Intervals(c(-Inf,lb_max)))
  }
  if(ub_min<Inf){
    exclude_intervals=interval_union(exclude_intervals,Intervals(c(ub_min, Inf)))
  }
  
  se=sqrt(sum(eta^2))*sqrt(sigma2)
  value=c(t(eta)%*%y)
  
  #to find CI, we need to find L,U such that F_(mean=L, sd=se)(eta^Ty)=alpha/2, F_(mean=U, sd=se)(eta^Ty)=1-alpha/2,
  ##where F_(a,b) is the cdf of Trunc t with mean a and sd b, with domain being complement of exclude_intervals
  f1=function(mu){return (cdfF(mu, se, df=n-sum(selected)-1, exclude_intervals,value)-alpha/2)}
  f2=function(mu){return (cdfF(mu, se, df=n-sum(selected)-1, exclude_intervals,value)-1+alpha/2)}
  L=nleqslv(value+qnorm(alpha/2)*se,f2)$x
  U=nleqslv(value+qnorm(1-alpha/2)*se, f1)$x
  return(c(L, U))
  
}

#cdf F(t) of truncated normal, given mean=mu, sd=se, and domain intervals
##property: F is monotone in its mean mu
cdfF=function(mu, se, df=0, exclude_intervals,value){
  #input:
  #- mu: mean parameter in truncNormal
  #- se: standard error parameter in truncNormal
  #- exclude_intervals: truncated-out region in truncNormal
  #- value: observed value of the truncNormal random variable
  
  #output:
  #- value of cdf F(value) of truncated normal
  
  normalizer=1-sum(pnorm(exclude_intervals@.Data[,2], mean = mu, sd=se)-pnorm(exclude_intervals@.Data[,1], mean=mu, sd=se))
  int_left_endpoints=pnorm(exclude_intervals@.Data[,1], mean = mu,sd=se)
  if(length(int_left_endpoints)>1){
    int_left_endpoints=int_left_endpoints-c(0, pnorm(exclude_intervals@.Data[-length(int_left_endpoints),2], mean = mu, sd=se))
    int_left_endpoints=cumsum(int_left_endpoints)
  }
  int_left_endpoints=int_left_endpoints/normalizer
  value_loc=which(exclude_intervals@.Data[,1]>value)[1]-1
  if(is.na(value_loc)){value_loc=length(int_left_endpoints)}
  result=int_left_endpoints[value_loc]+(pnorm(value, mean = mu, sd=se)-pnorm(exclude_intervals@.Data[value_loc,2], mean = mu, sd=se))/normalizer
  
  return(result)
}


#helper func: find the quadratic root
quad <- function(A,B,C){
  #input: in quadratic functions A*x^2 + B*x + C = 0
  #- A: coefficient of the quadratic term x^2
  #- B: coefficient of the linear term x
  #- C: the constant
  
  #output:
  #- the solution to this quadratic function
  
  answer <- c((-B - sqrt(B^2 - 4 * A * C)) / (2 * A),
              (-B + sqrt(B^2 - 4 * A * C)) / (2 * A))
  return(answer)
}

# Get standard error of mean response for new data samples, using ground truth sigma=1
# An element at index i is a standard error of the ith data sample.
# Parameters: X is the data used to fit the model; new_x is the new data points
get_standard_errors <- function(X,new_x) {
  #input:
  #- X: design matrix
  #- new_x: a new x point
  
  #output:
  #- standard error of the predicted mean at new_x
  
  cov_inv <- t(X) %*%  X
  cov_matrix <- solve(cov_inv)
  std_error_matrix <- new_x %*% cov_matrix %*% t(new_x)
  std_error=sqrt(diag(std_error_matrix))
  std_error
}

#helper func: to facilitate parallel computation for each dataset
#run model selection with the specified model selection criteria, and return 95% CI for response mean
subfunc=function(statistic='aic', all_compete_models=NA, alpha=0.05,sigma2){
  #input: 
  #- statistic: model selection criteria, must be one of aic, bic, aicc
  #- all_compete_models: binary representation of all submodels to be compared
  #- alpha: significance level
  #- sigma2: vairance parameter if known
  
  #output:
  #- classical CI for the predicted mean at new x points when sigma2 is known(based on normal distribution),
  #post-selection CI, the selected model size, and the selected model

  if(is.na(all_compete_models)){
    all_compete_models=(alls==1)
  }
  
  all_vars <- c(paste("X", 1:p, sep=""))
  
  df=data.frame(X)
  y <- X[,1:length(B)]%*%B + rnorm(nrow(X), 0, sqrt(sigma2))
  df$Y=y
  # Choose the best subset of variables.
  models <- regsubsets(Y~., data = df, nvmax = p, method = 'exhaustive', nbest = 1, intercept=TRUE)
  reg_summary = summary(models)
  
  if (statistic == 'aic') {
    aics <- rep(NA,ncol(X))
    for (i in 1:ncol(X)) {
      vars <- all_vars[reg_summary$which[i,-1]]
      # Fit the linear regression with the variables in subset.
      fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
      est_fit <- lm(fmla, data = df)
      aics[i] <- AIC(est_fit)
    }
    size = which.min(aics)
  } else if (statistic == 'bic') {
    size = which.min(reg_summary$bic)
  } else if(statistic=='aicc'){
    aiccs <- rep(NA,ncol(X))
    for (i in 1:ncol(X)) {
      vars <- all_vars[reg_summary$which[i,-1]]
      # Fit the linear regression with the variables in subset.
      fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
      est_fit <- lm(fmla, data = df)
      aiccs[i] <- AIC(est_fit)+2*length(vars)*(length(vars)+1)/(n-length(vars)-1)
    }
    size = which.min(aiccs)
  }else {
    print("Invalid statistic!")
    break;
  }
  
  # Create an array of variables in the subset
  vars <- all_vars[reg_summary$which[size,-1]]
  
  # Fit the linear regression with the variables in subset
  fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
  est_fit <- lm(fmla, data = df)
  #Get predisction as well as CIs(using the estimated variance) for response mean of selected models
  est_ci_newdt=predict(est_fit, data.frame(new_x), interval = "confidence",level=1-alpha)
  
  
  # Get CIs for response mean of selected model using the true variance.
  
  # Get a vector of standard errors for selected model.
  # Generate a matrix containing only selected variables.
  data_est <- rep(1,nrow(X)) # Create the first column for alignment.
  newdt_est=rep(1,num_newdt)
  for (i in 1:ncol(X)) {
    if (all_vars[i] %in% vars) {
      data_est <- cbind(data_est, X[,i])
      newdt_est=cbind(newdt_est,new_x[,i])
    }
  }
  std_error_est <- get_standard_errors(data_est,newdt_est)*sqrt(sigma2)
  est_lwr_var <- est_ci_newdt[,1] - qnorm(alpha/2,lower.tail = F) * std_error_est
  est_upr_var <- est_ci_newdt[,1] + qnorm(alpha/2,lower.tail = F) * std_error_est
  
  
  
  #get corrected CIs for response mean using our method
  lur_bounds=apply(X=new_x, MARGIN=1, FUN=postICci, X_dt=X, y=df$Y, selected=reg_summary$which[size,-1], 
                   alls=all_compete_models, criteria=statistic,alpha=alpha,sigma2=sigma2)
  est_lwr_corrected=lur_bounds[1,]
  est_upr_corrected=lur_bounds[2,]
  
  result=c(est_lwr_var, est_upr_var, est_lwr_corrected, est_upr_corrected, size, paste(vars,collapse = '+'))
  return(result)
}


# the simulation fucntion: parallel version------------
simulate <- function(N, B, p, statistic, sigma2=1,alpha=0.05){
  all_vars <- c(paste("X", 1:p, sep=""))
  overfitted = 0
  
  # #record the coverage of the IC-selecte model using our corrected CI
  # est_ci_corrected_coverage = rep(0,num_newdt)
  # #record the coverage of the IC-selecte model, with true/estimated sigma2, but no correction
  # est_ci_coverage = rep(0,num_newdt)
  
  
  #record the coverage freq for different sizes (from 1 to p) of selected model
  coverage_selected=matrix(0,nrow = p,ncol = num_newdt)
  #record the coverage freq for different sizes (from 1 to p) of selected model, after our correction
  coverage_selected_corrected=coverage_selected
  
  mean_ynew=new_x[,1:length(B)]%*%B
  #build a dataframe
  df=data.frame(X)
  
  #record cases when numerical issue happens(i.e. normalizer=0)
  numerical_error_count=rep(0,nrow(new_x))
  
  #facilitating variables in our correction 
  z=expand.grid(rep(list(0:1),p))[-1,]
  alls=matrix(unlist(z), ncol = p, byrow = F)
  
  #registerDoParallel(25)
  all_results=foreach(j=1:N,.combine = rbind)%dopar%{
    set.seed(j)
    subfunc(statistic=statistic, all_compete_models = (alls==1), alpha = alpha, sigma2 = sigma2)
  }
  
  num_newdt=nrow(new_x)
  mean_ynew=c(mean_ynew)
  #record the selected models
  mdls_selected=all_results[,ncol(all_results)]
  overfitted=sum(grepl('X1+X2+X3+X', mdls_selected, fixed = TRUE))
  
  #calculate the coverage
  all_results=matrix(as.numeric(all_results[,-ncol(all_results)]), nrow = nrow(all_results))
  est_ci_coverage=(t(t(all_results[,1:num_newdt])<=mean_ynew) & t(t(all_results[,(num_newdt+1):(num_newdt*2)])>=mean_ynew))
  est_ci_corrected_coverage=(t(t(all_results[,(num_newdt*2+1):(num_newdt*3)])<=mean_ynew) & t(t(all_results[,(num_newdt*3+1):(num_newdt*4)])>=mean_ynew))
  
  all_selected_sizes=as.integer(all_results[,ncol(all_results)])
  size_selected=colSums(sapply(1:p,FUN='==',all_selected_sizes))#count of selected model size
  
  
  for (i in 1:p) {
    rid=which(all_selected_sizes==i)
    if(length(rid)==1){
      coverage_selected[i,]=est_ci_coverage[rid,]*1
      coverage_selected_corrected[i,]=est_ci_corrected_coverage[rid,]*1
    }else if(length(rid)==0){
      coverage_selected[i,]=NA
      coverage_selected_corrected[i,]=NA
    }else{
      coverage_selected[i,]=colSums(est_ci_coverage[rid,])
      coverage_selected_corrected[i,]=colSums(est_ci_corrected_coverage[rid,])
    }
  }
  
  print(paste("number of strictly overfitted model:",overfitted))
  print("frequency of selected model size (from 1-p accordingly):")
  print(size_selected)
  
  print("Using classical methods, overall coverage at each new_xpoints:")
  print(colMeans(est_ci_coverage))
  
  print("coverage of selected models at all new_x points, group by selected model size:")
  print(coverage_selected/matrix(rep(size_selected,num_newdt),nrow = p,byrow = F))
  
  print("Using corrected methods, overall coverage at each new_xpoints:")
  print(colMeans(est_ci_corrected_coverage))
  
  print("After correction, coverage of selected models at all new_x points, group by selected model size:")
  print(coverage_selected_corrected/matrix(rep(size_selected,num_newdt),nrow = p,byrow = F))
  
}


# Function generates a correlation matrix corresponding to AR(1) dependence structure.
ar1_matrix <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n-1))
  rho^exponent
}

# Function generates an equal correlation matrix.
equal_cor_matrix <- function(n, rho) {
  rho*matrix(rep(1, n ^ 2), nrow = n, ncol = n) + (1-rho)*diag(n)
}


#===============PART2: Simulation settings==============================
N <- 5000 #num of iterations
n=50 #sample size
p=10 #num of candidate predictors
B<- c(1, 2, 3) #ground truth coefficient(beta)
sigma2=1
# Possible statistics: "aic", "aicc", "bic"
statistic <- "aic"

#-specify covariance matrix for the design matrix X
# Possible scenarios: 1, 2, 3, 4, 5
# scenario is an integer corresponding to the following:
# 1: AR(1) with rho = 0.5
# 2: AR(1) with rho = -0.5
# 3: IID error
# 4: equal correlation with rho = 0.5 
# 5: equal correlation with rho = -0.5 
scenario <- 1
cor_matrix <- switch(scenario, ar1_matrix(p, 0.5), ar1_matrix(p, -0.5), 
                     diag(p), equal_cor_matrix(p, 0.5), equal_cor_matrix(p, -0.5))


#===============PART3: EXPERIMENTS==============================
set.seed(1)
# some new x points to test performance of prediction
num_newdt=10 #number of new data points whose predictive performance is to be tested
new_x=mvrnorm(n=num_newdt, mu = rep(0, p), Sigma = cor_matrix)
set.seed(2)
X=mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)
set.seed(3)
start_time <- Sys.time()
simulate(N=N, B=B, p=p, statistic='aic',alpha = 0.05,sigma2 = sigma2)
end_time <- Sys.time()
print(end_time-start_time)

#===output:
# [1] "number of strictly overfitted model: 3662"

# [1] "frequency of selected model size (from 1-p accordingly):"
# [1]    0    0 1337 1537 1203  615  241   61    5    1

# [1] "Using classical methods, overall coverage at each new_xpoints:"
# [1] 0.9048 0.9014 0.8682 0.9328 0.8738 0.9070 0.9024 0.8526 0.9182 0.8744

# [1] "coverage of selected models at all new_x points, group by selected model size:"
# [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]     [,10]
# [1,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [2,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [3,] 0.9491399 0.9424084 0.9528796 0.9573672 0.9439043 0.9506358 0.9521316 0.9521316 0.9483919 0.9424084
# [4,] 0.9316851 0.9186727 0.8835394 0.9264802 0.9089135 0.9232271 0.9284320 0.8347430 0.9037085 0.9082628
# [5,] 0.8894431 0.8935993 0.8320865 0.9110557 0.8445553 0.8861180 0.8861180 0.8021613 0.9052369 0.8561929
# [6,] 0.8308943 0.8292683 0.7951220 0.9414634 0.7544715 0.8487805 0.8162602 0.8292683 0.9121951 0.7804878
# [7,] 0.8049793 0.8257261 0.7178423 0.9211618 0.7593361 0.8381743 0.8049793 0.7634855 0.9253112 0.7012448
# [8,] 0.6885246 0.7540984 0.7049180 0.9508197 0.6885246 0.8196721 0.7377049 0.7049180 0.9180328 0.5409836
# [9,] 1.0000000 0.8000000 0.6000000 1.0000000 0.8000000 1.0000000 0.8000000 0.8000000 0.8000000 0.6000000
# [10,] 1.0000000 1.0000000 0.0000000 0.0000000 1.0000000 0.0000000 1.0000000 1.0000000 1.0000000 1.0000000

# [1] "Using corrected methods, overall coverage at each new_xpoints:"
# [1] 0.9490 0.9454 0.9486 0.9468 0.9442 0.9470 0.9468 0.9440 0.9486 0.9466

# [1] "After correction, coverage of selected models at all new_x points, group by selected model size:"
# [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]     [,10]
# [1,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [2,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [3,] 0.9491399 0.9424084 0.9528796 0.9573672 0.9439043 0.9506358 0.9521316 0.9521316 0.9483919 0.9424084
# [4,] 0.9577098 0.9479506 0.9577098 0.9466493 0.9472999 0.9420950 0.9505530 0.9427456 0.9472999 0.9603123
# [5,] 0.9426434 0.9443059 0.9393184 0.9260183 0.9376559 0.9501247 0.9443059 0.9434746 0.9517872 0.9459684
# [6,] 0.9414634 0.9577236 0.9560976 0.9626016 0.9495935 0.9463415 0.9398374 0.9479675 0.9495935 0.9349593
# [7,] 0.9377593 0.9294606 0.9045643 0.9460581 0.9419087 0.9502075 0.9336100 0.9045643 0.9460581 0.9211618
# [8,] 0.9672131 0.9016393 0.9180328 0.9836066 0.9508197 0.9344262 0.9016393 0.9180328 0.9344262 0.9180328
# [9,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 0.8000000 1.0000000
# [10,] 1.0000000 1.0000000 0.0000000 0.0000000 1.0000000 0.0000000 1.0000000 1.0000000 1.0000000 1.0000000
