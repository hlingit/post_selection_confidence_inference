#replicate the simulation in "On overfitting and post-selection uncertainty assessments"
#for motivation example in our paper
rm(list=ls())
library(foreach)
library(leaps)
library(MASS)
source('sourceCode.R')

#========PART 1: helper func for data generation and model selection========
#check if the model is strictly overfitted
overfit <- function(vars) {
  #input: a sequence of selected variable names (in the form of 'X1','X2',...)
  
  #output:assuming the oracle model contains only 'X1','X2','X3',
  #- a boolean value, indicating whether the input contains more variables than 'X1','X2','X3'

  true_vars <- c("X1", "X2", "X3")
  all(true_vars %in% vars) & (length(vars) > 3)
}

#run model selection with the specified model selection criteria, and return 95% CI for response mean
subfunc=function(statistic='aic', all_compete_models=NA, alpha=0.05,sigma2){
  #input: 
  #- statistic: model selection criteria, must be one of aic, bic, aicc
  #- all_compete_models: binary representation of all submodels to be compared
  #- alpha: significance level
  #- sigma2: vairance parameter if known
  
  #output:
  #- CI for the predicted mean at new x points when sigma2 is unknown (based on t-distribution),
  #when sigma2 is known(based on normal distribution), the selected model size, and the selected model
  
  all_vars <- c(paste("X", 1:p, sep=""))
  
  df=data.frame(X)
  y <- X[,1:length(B)]%*%B + rnorm(nrow(X), 0, 1)
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
  } else{
    print("Invalid statistic!")
    break;
  }
  
  
  # Create an array of variables in the subset
  vars <- all_vars[reg_summary$which[size,-1]]
  
  # Fit the linear regression with the variables in subset
  fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
  est_fit <- lm(fmla, data = df)
  #Get predisction as well as CIs(using the estimated variance) for response mean of selected models
  est_ci_newdt=predict(est_fit, data.frame(new_x), interval = "confidence")
  
  # Get a vector of standard errors for selected model.
  # Generate a matrix containing only selected variables.
  data_est <- rep(1,nrow(X)) # Create the first column for alignment.
  newdt_est=rep(1,num_newdt)
  for (i in 1:p) {
    if (all_vars[i] %in% vars) {
      data_est <- cbind(data_est, X[,i])
      newdt_est=cbind(newdt_est,new_x[,i])
    }
  }
  #data_est <- data_est[,-1] # Delete the first column for alignment.
  std_error_est <- get_standard_errors(data_est,newdt_est)
  
  # Get CIs for response mean of selected model using the true variance.
  est_lwr_truevar <- est_ci_newdt[,1] - qnorm(alpha/2,lower.tail = F) * std_error_est
  est_upr_truevar <- est_ci_newdt[,1] + qnorm(alpha/2,lower.tail = F) * std_error_est
  
  result=c(est_ci_newdt[,2], est_ci_newdt[,3], est_lwr_truevar, est_upr_truevar, size, paste(vars,collapse = '+'))
  return(result)
}

#simulation function
simulate <- function(N, B, p, statistic, alpha=0.05,sigma2=1){
  all_vars <- c(paste("X", 1:p, sep=""))
  overfitted = 0
  
  #record the coverage of the IC-selecte model using classical t distribution-based  CI
  est_ci_t_coverage = rep(0,num_newdt)
  #record the coverage of the IC-selecte model, with true sigma2, but no correction
  est_ci_true_coverage = rep(0,num_newdt)
  #record the selected size
  size_selected=rep(0,p)
  
  #record the coverage freq for different sizes (from 1 to p) of selected model, using sigma=1
  coverage_true_selected=matrix(0,nrow = p,ncol = num_newdt)
  #record the coverage freq for different sizes (from 1 to p) of selected model, with t distribution-based CI
  coverage_t_selected=coverage_true_selected
  
  mean_ynew=new_x[,1:length(B)]%*%B
  
  #matrix represenation of all possible models from candidate predictors
  z=expand.grid(rep(list(0:1),p))[-1,]
  alls=matrix(unlist(z), ncol = p, byrow = F)
  
  #registerDoParallel(25)
  all_results=foreach(j=1:N,.combine = rbind)%dopar%{
    set.seed(j)
    subfunc(statistic=statistic, all_compete_models = (alls==1), alpha = alpha, sigma2 = sigma2)
  }
  
  
  #record the selected models
  mdls_selected=all_results[,ncol(all_results)]
  overfitted=sum(grepl('X1+X2+X3+X', mdls_selected, fixed = TRUE))
  #calculate the coverage
  all_results=matrix(as.numeric(all_results[,-ncol(all_results)]), nrow = nrow(all_results))
  mean_ynew=c(mean_ynew)
  est_ci_t_coverage=(t(t(all_results[,1:num_newdt])<=mean_ynew) & t(t(all_results[,(num_newdt+1):(num_newdt*2)])>=mean_ynew))
  est_ci_true_coverage=(t(t(all_results[,(num_newdt*2+1):(num_newdt*3)])<=mean_ynew) & t(t(all_results[,(num_newdt*3+1):(num_newdt*4)])>=mean_ynew))
  
  all_selected_sizes=as.integer(all_results[,ncol(all_results)])
  size_selected=colSums(sapply(1:p,FUN='==',all_selected_sizes))#count of selected model size
  
  
  print(paste("number of strictly overfitted model:",overfitted))
  print("frequency of selected model size (from 1-p accordingly):")
  print(size_selected)
  
  print("Using sigma=1, overall coverage at each new_xpoints:")
  print(colMeans(est_ci_true_coverage))
  
  print("With unkown sigma, overall coverage at each new_xpoints:")
  print(colMeans(est_ci_t_coverage))
  
}

#==========PART 2: model parameter setting===========
#num of predictors
p=10
#sample size
n=50
#ground truth beta
B <- c(1, 2, 3)
#model selection criteria
statistic <- "aic"
#1st-order auto-regressive covariance matrix with correlation parameter
cor_matrix=ar1_matrix(p, 0.5)
#num of new x points whose 95% CI will be constructed for its predicted mean
num_newdt=10
set.seed(1)
# generate the new x
new_x=MASS::mvrnorm(n=num_newdt, mu = rep(0, p), Sigma = cor_matrix)
set.seed(2)
# generate the design matrix
X=MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)
set.seed(3)


#===========PART 3: run simulation==========
simulate(N=5000,B,p,statistic = 'aic')

# [1] "number of strictly overfitted model: 3662"
# [1] "frequency of selected model size (from 1-p accordingly):"
# [1]    0    0 1337 1537 1203  615  241   61    5    1
# [1] "Using sigma=1, overall coverage at each new_xpoints:"
# [1] 0.9048 0.9014 0.8682 0.9328 0.8738 0.9070 0.9024 0.8526 0.9182 0.8744
# [1] "With unkown sigma, overall coverage at each new_xpoints:"
# [1] 0.8916 0.8896 0.8572 0.9268 0.8636 0.8964 0.8920 0.8394 0.9088 0.8602


