#a textbook example from https://otexts.com/fpp3/regression-intro.html, chap 7.5

rm(list = ls())
source('/Users/linhm/Desktop/Li/postAIC/code for github/sourceCode.R')
# load packages
library(leaps)
library(intervals)
library(natural)


#==============PART 1: Helper functions===============
# perform post-selection correction for confidence inference
postICci_general=function(X_dt, y, selected, alls,eta, criteria='aic', alpha=0.05){
  #input:
  #- X_dt: design matrix
  #- y: response
  #- selected: the AIC selected model (boolean vector indicating whether the predictor is included or not)
  #- alls: all other models to be compard with the selected model (binary matrix: each row represent a model)
  #- criteria: model selection criteria, must be one of 'aic', 'bic', 'aicc'
  #- alpha: significance level
  #- eta: a vector specified by the user; eta^TY is the research target. Some examples of eta:
  #--when the target is prediction mean at a new data point, then eta=cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%c(1,new_xpoint[selected])
  #--when the target is each coefficient: then eta is each row of solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
  
  n=length(y)
  if(criteria=='aic'){
    cons=2
  }else if(criteria=='bic'){
    cons=log(n)
  }else if(criteria=='aicc'){
    # pass
  }else{
    print("must be either aic or bic!")
    return(NA)
  }
  #record the excluded intervals for those TN with 2 intervals
  exclude_lb=rep(NA,nrow(alls))
  exclude_ub=exclude_lb
  #record the included intervals for those TN with 1 interval
  lb_max=-Inf
  ub_min=Inf
  
  P_select=diag(n)-cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
  cee=eta/c(t(eta)%*%eta)
  z=y-cee%*%t(eta)%*%y
  
  for (s in 1:nrow(alls)) {
    #skip comparing with the selected model itself
    if(sum(alls[s,]!=selected)==0){
      next
    }
    
    #exp weight involving cons
    s_hat=sum(selected)#size of selected model
    s_tp=sum(alls[s,])#size of model alls[s,]
    if(criteria=='aic' | criteria=='bic'){
      wt=exp(cons*(s_hat-s_tp)/n)
    }
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
  
  #if sigmahat is not provided, then we need to estimate
  #mse-based estiamte
  sigmahat1=sigma(lm(y~X_dt))
  #organic lasso estiamte
  sigmahat2=olasso_cv(x = X_dt, y =y)$sig_obj
  
  se1=sqrt(sum(eta^2))*sigmahat1
  se2=sqrt(sum(eta^2))*sigmahat2
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
  return(c(L1,U1,L2,U2))
  
}

#=============PART 2: Read in data===============
us_change=read.csv('/Users/linhm/Desktop/Li/postAIC/us_change.csv')
us_change=us_change[,2:6]

#num of candidate predictors
p=4
#sample size
n=nrow(us_change)
#names of the p candidate predictors
all_vars=c('Income','Production','Savings','Unemployment')
models <- regsubsets(Consumption~Income+Production+Savings+Unemployment, data = us_change, 
                     nvmax =p, method = 'exhaustive', nbest = 1, intercept=TRUE)
reg_summary = summary(models)

# find the aic selected model
aics <- rep(NA,p)
for (i in 1:p) {
  vars <- all_vars[reg_summary$which[i,-1]]
  # Fit the linear regression with the variables in subset.
  fmla <- as.formula(paste("Consumption~", paste(vars, collapse="+")))
  est_fit <- lm(fmla, data = us_change)
  aics[i] <- AIC(est_fit)
}
size = which.min(aics)
# aic selected variables
vars <- all_vars[reg_summary$which[size,-1]]
mdl_aic=lm(as.formula(paste("Consumption~", paste(vars, collapse="+"))), data = us_change)
# classical 95% CIs for each coef
tp=summary(mdl_aic)$coef[-1,]
cis_naive=cbind(tp[,1]+qt(0.05/2, df=n-p-1)*tp[,2], tp[,1]+qt(1-0.05/2, df=n-p-1)*tp[,2])
colnames(cis_naive)=c('l0','u0')
selected=reg_summary$which[size,-1]
print('Classical 95% CI for each regression coefficient:')
cis_naive

#===========PART 3: 95% CIs with post-aic correction=============
# convert the df into matrix
X_dt=data.matrix(us_change[,-1], rownames.force = NA)
# eta of the point estimator eta^TY 
etas=solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
#all models that are compared during the aic procedure
z=expand.grid(rep(list(0:1),p))[-1,]
alls=matrix(unlist(z), ncol = p, byrow = F)
alls=alls==1
cis_corrected=matrix(nrow = sum(selected),ncol = 4)
for (i in 1:sum(selected)) {
  cis_corrected[i,]=postICci_general(X_dt=X_dt,y=us_change$Consumption,selected=selected,alls=alls,eta=etas[i+1,], criteria='aic',alpha=0.05)
  
}
colnames(cis_corrected)=c('l1_mse','u1_mse','l2_olasso','u2_olasso')
rownames(cis_corrected)=rownames(cis_naive)
print('Post-AIC corrected 95% CI for each regression coefficient:')
cis_corrected


#===output:
# [1] "Classical 95% CI for each regression coefficient:"
# > cis_naive
#                        l0          u0
# Income        0.661463330  0.81970365
# Production    0.001528875  0.09281636
# Savings      -0.058657459 -0.04712279
# Unemployment -0.363064140  0.01369361
# 
# [1] "Post-AIC corrected 95% CI for each regression coefficient:"
# > cis_corrected
#                   l1_mse      u1_mse   l2_olasso   u2_olasso
# Income        0.66146524  0.82526995  0.65446730  0.83546913
# Production   -0.01115300  0.11520277 -0.01378375  0.11756396
# Savings      -0.05938233 -0.04712340 -0.06023854 -0.04661389
# Unemployment -0.45825078  0.06364616 -0.46776089  0.07568029