# To replicate the simu results in sec 4 of the paper, with unknown sigma2.

#Overview: Three strategies for estimating sigma
#---strategy1: mse based on full model
#---strategy2: organic lasso
#---strategy3: mse based on aic-selected model

rm(list = ls())
#load packages
library(intervals)#library(sets)
library(nleqslv)
library(leaps)
library(natural)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)
library(bigstatsr)

#===================PART1: helper functions============================
# perform post-selection correction for confidence inference
postICci=function(X_dt, y, selected, alls, new_xpoint, statistic='aic', alpha=0.05){
  #input:
  #- X_dt: design matrix
  #- y: response
  #- selected: the selected model (vector: boolean element indicating whether the predictor is included or not)
  #- alls: all other models to be compard with the selected model (binary matrix: each row represent a model)
  #- new_xpoint: new x point
  #- statistic: model selection criteria, must be one of 'aic', 'bic', 'aicc'
  #- alpha: significance level
  #- sigma2: ground truth noise level in the model
  
  #output:
  #- the post-selection confidence interval for predicted mean at new_xpoint
  
  
  n=length(y)
  if(statistic=='aic'){
    cons=2
  }else if(statistic=='bic'){
    cons=log(n)
  }else if(statistic=='aicc'){
    # pass
  }
  else{
    print("must be either aic, aicc or bic!")
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
    #skip comparing with the selected model itself or its superset
    if(all(alls[s,]-selected>=0)){
      next
    }
    
    #exp weight involving cons
    s_hat=sum(selected)#size of selected model
    s_tp=sum(alls[s,])#size of model alls[s,]
    if(criteria=='aic' | criteria=='bic'){
      wt=exp(cons*(s_hat-s_tp)/n)
    }
    if(statistic=='aicc'){
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
  
  #if sigmahat is not provided, then we need to estimate
  #mse-based estiamte from the full model
  sigmahat1=sigma(lm(y~X_dt))
  #organic lasso estiamte
  sigmahat2=olasso_cv(x = X_dt, y =y)$sig_obj
  #mse-based estimate from the selected model
  sigmahat3=sigma(lm(y~X_dt[,selected]))
  
  se1=sqrt(sum(eta^2))*sigmahat1
  se2=sqrt(sum(eta^2))*sigmahat2
  se3=sqrt(sum(eta^2))*sigmahat3
  value=c(t(eta)%*%y)
  
  #to find CI, we need to find L,U such that F_(mean=L, sd=se)(eta^Ty)=alpha/2, F_(mean=U, sd=se)(eta^Ty)=1-alpha/2,
  ##where F_(a,b) is the cdf of Trunc t with mean a and sd b, with domain being complement of exclude_intervals
  f1=function(mu){return (cdfF(mu, se1, df=n-sum(selected)-1, exclude_intervals,value)-alpha/2)}
  f2=function(mu){return (cdfF(mu, se1, df=n-sum(selected)-1, exclude_intervals,value)-1+alpha/2)}
  L1=nleqslv(value+qnorm(alpha/2)*se1,f2)$x
  U1=nleqslv(value+qnorm(1-alpha/2)*se1, f1)$x
  f1=function(mu){return (cdfF(mu, se2, df=n-sum(selected)-1, exclude_intervals,value)-alpha/2)}
  f2=function(mu){return (cdfF(mu, se2, df=n-sum(selected)-1, exclude_intervals,value)-1+alpha/2)}
  L2=nleqslv(value+qnorm(alpha/2)*se2,f2)$x
  U2=nleqslv(value+qnorm(1-alpha/2)*se2, f1)$x
  f1=function(mu){return (cdfF(mu, se3, df=n-sum(selected)-1, exclude_intervals,value)-alpha/2)}
  f2=function(mu){return (cdfF(mu, se3, df=n-sum(selected)-1, exclude_intervals,value)-1+alpha/2)}
  L3=nleqslv(value+qnorm(alpha/2)*se3,f2)$x
  U3=nleqslv(value+qnorm(1-alpha/2)*se3, f1)$x
  return(c(L1,U1,L2,U2,L3,U3))
  
}

# cdf F(t) of Truncated Normal, given mean=mu, sd=se, and domain
cdfF=function(mu, se, df=0, exclude_intervals, value){
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


# helper func: find the quadratic root
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



#helper func: run model selection with the specified model selection criteria, and return 95% CI for response mean
subfunc=function(statistic='aic', all_compete_models=NA, alpha=0.05, new_x){
  #input: 
  #- statistic: model selection criteria, must be one of aic, bic, aicc
  #- all_compete_models: binary matrix representation of all submodels to be compared
  #- alpha: significance level
  #- new_x: the new x data points at which the predicted mean needs a confidence interval
  
  #output:
  #- post-selection CIs for the predicted mean at new x points when sigma2 is unknown using the 3 
  #estimating strategies for sigma, the selected model size, and the selected model
  
  all_vars <- c(paste("X", 1:p, sep=""))
  
  df=data.frame(X)
  y <- X[,1:length(B)]%*%B + rnorm(nrow(X), 0, sqrt(1))#sigma2=1 in data generation
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
  
  
  #get corrected CIs for response mean using our method
  lur_bounds=apply(X=new_x, MARGIN=1, FUN=postICci, X_dt=X, y=df$Y, selected=reg_summary$which[size,-1], 
                   alls=all_compete_models, statistic=statistic,alpha=alpha)
  #CI from mse-based sigmahat from full model
  est_lwr_corrected1=lur_bounds[1,]
  est_upr_corrected1=lur_bounds[2,]
  #CI from organic lasso sigmahat
  est_lwr_corrected2=lur_bounds[3,]
  est_upr_corrected2=lur_bounds[4,]
  #CI from mse-based sigmahat from selected model
  est_lwr_corrected3=lur_bounds[5,]
  est_upr_corrected3=lur_bounds[6,]
  
  result=c(est_lwr_corrected1, est_upr_corrected1,est_lwr_corrected2, est_upr_corrected2, 
           est_lwr_corrected3, est_upr_corrected3, size, paste(vars,collapse = '+'))
  return(result)
}


# simulation function: parallel version------------
simulate <- function(X, new_x, N, B, p,statistic, alpha=0.05){
  all_vars <- c(paste("X", 1:p, sep=""))
  overfitted = 0
  
  #record the coverage freq for different sizes (from 1 to p) of selected model, after our correction
  coverage_selected_corrected1=matrix(0,nrow = p,ncol = num_newdt)
  coverage_selected_corrected2=matrix(0,nrow = p,ncol = num_newdt)
  coverage_selected_corrected3=matrix(0,nrow = p,ncol = num_newdt)
  
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
    subfunc(statistic=statistic, all_compete_models = (alls==1), alpha = alpha, new_x = new_x)
    #subfunc(statistic='bic', all_compete_models = (alls==1),alpha = 0.05,sigmaKnown = F)
  }
  
  num_newdt=nrow(new_x)
  mean_ynew=c(mean_ynew)
  #record the selected models
  mdls_selected=all_results[,ncol(all_results)]
  overfitted=sum(grepl('X1+X2+X3+X', mdls_selected, fixed = TRUE))
  
  #calculate the coverage
  all_results=matrix(as.numeric(all_results[,-ncol(all_results)]), nrow = nrow(all_results))
  #coverage of CI from mse-based sigmahat from full model
  est_ci_corrected1_coverage=(t(t(all_results[,1:num_newdt])<=mean_ynew) & t(t(all_results[,(num_newdt+1):(num_newdt*2)])>=mean_ynew))
  #coverage of CI from organic lasso sigmahat
  est_ci_corrected2_coverage=(t(t(all_results[,(num_newdt*2+1):(num_newdt*3)])<=mean_ynew) & t(t(all_results[,(num_newdt*3+1):(num_newdt*4)])>=mean_ynew))
  #coverage of CI from mse-based sigmahat from selected model
  est_ci_corrected3_coverage=(t(t(all_results[,(num_newdt*4+1):(num_newdt*5)])<=mean_ynew) & t(t(all_results[,(num_newdt*5+1):(num_newdt*6)])>=mean_ynew))
  
  all_selected_sizes=as.integer(all_results[,ncol(all_results)])
  size_selected=colSums(sapply(1:p,FUN='==',all_selected_sizes))#count of selected model size
  
  
  for (i in 1:p) {
    rid=which(all_selected_sizes==i)
    if(length(rid)==1){
      #coverage_selected[i,]=est_ci_coverage[rid,]*1
      coverage_selected_corrected1[i,]=est_ci_corrected1_coverage[rid,]*1
      coverage_selected_corrected2[i,]=est_ci_corrected2_coverage[rid,]*1
      coverage_selected_corrected3[i,]=est_ci_corrected3_coverage[rid,]*1
    }else if(length(rid)==0){
      #coverage_selected[i,]=NA
      coverage_selected_corrected1[i,]=NA
      coverage_selected_corrected2[i,]=NA
      coverage_selected_corrected3[i,]=NA
    }else{
      #coverage_selected[i,]=colSums(est_ci_coverage[rid,])
      coverage_selected_corrected1[i,]=colSums(est_ci_corrected1_coverage[rid,])
      coverage_selected_corrected2[i,]=colSums(est_ci_corrected2_coverage[rid,])
      coverage_selected_corrected3[i,]=colSums(est_ci_corrected3_coverage[rid,])
    }
  }
  
  print(paste("number of strictly overfitted model:",overfitted))
  print("frequency of selected model size (from 1-p accordingly):")
  print(size_selected)
  
  print("Using corrected methods and full-model mse sigmahat, overall coverage at each new_xpoints:")
  print(colMeans(est_ci_corrected1_coverage))
  print("Using corrected methods and olasso sigmahat, overall coverage at each new_xpoints:")
  print(colMeans(est_ci_corrected2_coverage))
  print("Using corrected methods and selected-model mse sigmahat, overall coverage at each new_xpoints:")
  print(colMeans(est_ci_corrected3_coverage))
  
  
  print("After correction and using full-model mse sigmahat, coverage of selected models at all new_x points, group by selected model size:")
  print(coverage_selected_corrected1/matrix(rep(size_selected,num_newdt),nrow = p,byrow = F))
  
  print("After correction and using olasso sigmahat, coverage of selected models at all new_x points, group by selected model size:")
  print(coverage_selected_corrected2/matrix(rep(size_selected,num_newdt),nrow = p,byrow = F))
  
  print("After correction and using selected-model mse sigmahat, coverage of selected models at all new_x points, group by selected model size:")
  print(coverage_selected_corrected3/matrix(rep(size_selected,num_newdt),nrow = p,byrow = F))
  
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


##===============PART2: Simulation settings==============================
N <- 5000 #num of iterations
n=50 #sample size
p=10
B<- c(1, 2, 3)

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


#===============PART3: Run experiments==============================
set.seed(1)
# some new x points to test performance of prediction
num_newdt=10 #number of new data points whose predictive performance is to be tested
new_x=mvrnorm(n=num_newdt, mu = rep(0, p), Sigma = cor_matrix)
set.seed(2)
X=mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)
set.seed(3)
start_time <- Sys.time()
simulate(X=X, new_x = new_x, N=N, B=B, p=p, statistic='aic',alpha = 0.05)
end_time <- Sys.time()
print(end_time-start_time)


#===output:
# [1] "number of strictly overfitted model: 3662"

# [1] "frequency of selected model size (from 1-p accordingly):"
# [1]    0    0 1337 1537 1203  615  241   61    5    1

# [1] "Using corrected methods and full-model mse sigmahat, overall coverage at each new_xpoints:"
# [1] 0.9446 0.9414 0.9468 0.9446 0.9406 0.9426 0.9442 0.9432 0.9456 0.9426

# [1] "Using corrected methods and olasso sigmahat, overall coverage at each new_xpoints:"
# [1] 0.9750 0.9708 0.9744 0.9740 0.9734 0.9738 0.9744 0.9714 0.9766 0.9756

# [1] "Using corrected methods and selected-model mse sigmahat, overall coverage at each new_xpoints:"
# [1] 0.9384 0.9328 0.9402 0.9364 0.9334 0.9356 0.9362 0.9360 0.9388 0.9358

# [1] "After correction and using full-model mse sigmahat, coverage of selected models at all new_x points, group by selected model size:"
# [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]     [,10]
# [1,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [2,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [3,] 0.9536275 0.9446522 0.9573672 0.9566193 0.9476440 0.9536275 0.9633508 0.9573672 0.9536275 0.9498878
# [4,] 0.9551074 0.9453481 0.9577098 0.9492518 0.9401431 0.9381913 0.9466493 0.9446975 0.9466493 0.9538061
# [5,] 0.9343308 0.9426434 0.9301746 0.9260183 0.9326683 0.9443059 0.9368246 0.9351621 0.9426434 0.9376559
# [6,] 0.9284553 0.9447154 0.9512195 0.9560976 0.9430894 0.9300813 0.9268293 0.9463415 0.9463415 0.9300813
# [7,] 0.9170124 0.9004149 0.9087137 0.9045643 0.9377593 0.9460581 0.9211618 0.8962656 0.9253112 0.9045643
# [8,] 0.9672131 0.8852459 0.8852459 0.9836066 0.9508197 0.9016393 0.8852459 0.9180328 0.9180328 0.8852459
# [9,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 0.6000000 1.0000000
# [10,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000

# [1] "After correction and using olasso sigmahat, coverage of selected models at all new_x points, group by selected model size:"
# [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]     [,10]
# [1,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [2,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [3,] 0.9932685 0.9895288 0.9955123 0.9917726 0.9902767 0.9947644 0.9940165 0.9932685 0.9917726 0.9895288
# [4,] 0.9856864 0.9720234 0.9778790 0.9739753 0.9804815 0.9713728 0.9746259 0.9772284 0.9746259 0.9876383
# [5,] 0.9584372 0.9584372 0.9542810 0.9551122 0.9542810 0.9667498 0.9625935 0.9551122 0.9659185 0.9617623
# [6,] 0.9512195 0.9626016 0.9691057 0.9739837 0.9626016 0.9560976 0.9658537 0.9577236 0.9756098 0.9577236
# [7,] 0.9543568 0.9502075 0.9543568 0.9668050 0.9585062 0.9626556 0.9543568 0.9336100 0.9626556 0.9502075
# [8,] 0.9508197 0.9344262 0.9672131 1.0000000 0.9672131 0.9508197 0.9344262 0.9508197 0.9672131 0.9180328
# [9,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
# [10,] 1.0000000 1.0000000 0.0000000 0.0000000 1.0000000 0.0000000 1.0000000 1.0000000 1.0000000 1.0000000

# [1] "After correction and using selected-model mse sigmahat, coverage of selected models at all new_x points, group by selected model size:"
# [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]      [,8]      [,9]     [,10]
# [1,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [2,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
# [3,] 0.9461481 0.9356769 0.9513837 0.9476440 0.9371728 0.9446522 0.9513837 0.9498878 0.9476440 0.9424084
# [4,] 0.9486012 0.9323357 0.9525049 0.9394925 0.9316851 0.9336370 0.9375407 0.9355888 0.9368900 0.9440468
# [5,] 0.9310058 0.9368246 0.9218620 0.9185370 0.9276808 0.9301746 0.9301746 0.9285121 0.9368246 0.9334996
# [6,] 0.9219512 0.9382114 0.9430894 0.9479675 0.9382114 0.9317073 0.9235772 0.9414634 0.9414634 0.9268293
# [7,] 0.9045643 0.9004149 0.9004149 0.9045643 0.9377593 0.9460581 0.9211618 0.8921162 0.9253112 0.8962656
# [8,] 0.9672131 0.8852459 0.8852459 0.9836066 0.9508197 0.9016393 0.8852459 0.9180328 0.9016393 0.8852459
# [9,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 0.6000000 1.0000000
# [10,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
