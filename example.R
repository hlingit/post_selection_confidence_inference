#confidence interval for predicted mean response, with post-AIC correction
#the noise level sigma can either be known or unknown
rm(list = ls())
source('sourceCode.R')

#===Basic setting===============
n=50 #sample size
p=10 # num of columns in X
B=rep(0,p)
B[1:3]=1:3 #ground truth model: only the first 3 predictors have nonzero effect on y
statistic <- "aic" #must be one of c("aic", "aicc", "bic")
alpha=0.05
sigmaKnown=F#assume noise level sigma is unknown


#===Data generation==============
# correlation structure of X: AR(1) with rho = 0.5
ar1_matrix <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n-1))
  rho^exponent
}
cor_matrix=ar1_matrix(p, 0.5)

#generate data used in model fitting
set.seed(1)
X=MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)
df=data.frame(X)
y=X%*%B + rnorm(nrow(X), 0, 1)
df$Y=y
#generate a new data point
new_x=MASS::mvrnorm(1, mu = rep(0, p), Sigma = cor_matrix)


#===Go through the information-based model selection procedure========
all_vars <- c(paste("X", 1:p, sep=""))
names(new_x)=all_vars
models <- leaps::regsubsets(Y~., data = df, nvmax = p, method = 'exhaustive', nbest = 1, intercept=TRUE)
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
#selected model
vars <- all_vars[reg_summary$which[size,-1]]#selected variables
fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
est_fit <- lm(fmla, data = df)
#estimated noise level, if unknown
if(!sigmaKnown){
  sigmahat=summary(est_fit)$sigma
}

#===Get post-selection corrected CI for predicted mean at new_x
#all possible models (use 1 to indicate a variable is selected and 0 otherwise)
z=expand.grid(rep(list(0:1),p))[-1,]
all_compete_models=(matrix(unlist(z), ncol = p, byrow = F))==1

corrected_CI=postICci(X, y, selected=reg_summary$which[size,-1], all_compete_models, new_x,
         criteria=statistic, alpha, sigmaKnown,sigmahat)
print(paste(1-alpha,' CI for predicted mean at new_x, with post-selection correction: (',
            round(corrected_CI[1],4),',',round(corrected_CI[2],4),')',sep = ''))
# output:
# "0.95 CI for predicted mean at new_x, with post-selection correction: (-7.4372,-5.7213)"

#In comparsion: before correction
est_ci_newdt=predict(est_fit, data.frame(t(new_x)), interval = "confidence",level=1-alpha)
print(paste(1-alpha,' CI for predicted mean at new_x before correction: (',
            round(est_ci_newdt[2],4),',',round(est_ci_newdt[3],4),')',sep = ''))
# output:
# "0.95 CI for predicted mean at new_x before correction: (-7.3661,-5.7212)"

#notice the 0.95 CI becomes wider (1.7159 vs 1.6449) after post-selection inference

