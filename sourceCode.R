library(intervals)
library(nleqslv)

postICci=function(X_dt,y,selected, alls, new_xpoint, criteria='aic',alpha=0.05,sigmahat=NA){
  #  Input Parameters:
  #  selected: the AIC selected model (vector: boolean element indicating whether the predictor is included or not)
  #  alls: all other models to be compard with AIC-selected model (matrix: each row represent a model)
  #  sigmahat: the ground truth noise level if sigma is Known=T, or an estimated value if sigma is unknown.
  
  #  Output:
  #  (1-alpha)-level confidence interval: c(L, U)
  
  n=length(y)
  if(criteria=='aic'){
    cons=2
  }else if(criteria=='bic'){
    cons=log(n)
  }else if(criteria=='aicc'){
    #pass; note: aicc=aic+(2k^2+2k)/(n-k-1)
  }else{
    print("must be either aic, aicc or bic!")
    return(NA)
  }
  
  #find the excluded intervals for those Truncated distributions
  exclude_lb=rep(NA,nrow(alls))
  exclude_ub=exclude_lb
  lb_max=-Inf
  ub_min=Inf
  
  eta=cbind(1,X_dt[,selected])%*%
    solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%
    c(1,new_xpoint[selected])
  cee=eta/c(t(eta)%*%eta)
  z=y-cee%*%t(eta)%*%y
  P_select=diag(n)-cbind(1,X_dt[,selected])%*%
    solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%
    t(cbind(1,X_dt[,selected]))
  
  for (s in 1:nrow(alls)) {
    #skip comparing with the selected model itself
    if(sum(alls[s,]!=selected)==0){
      next
    }
    s_hat=sum(selected)#size of selected model
    s_tp=sum(alls[s,])#size of model alls[s,]
    if(criteria=='aic' | criteria=='bic'){
      wt=exp(cons*(s_hat-s_tp)/n)
    }
    if(criteria=='aicc'){
      wt=exp((2*(s_hat-s_tp)+2*s_hat*(s_hat+1)/(n-s_hat-1)-2*s_tp*(s_tp+1)/(n-s_tp-1))/n)
    }
    P_s=diag(n)-cbind(1,X_dt[,alls[s,]])%*%
      solve(t(cbind(1,X_dt[,alls[s,]]))%*%cbind(1,X_dt[,alls[s,]]))%*%
      t(cbind(1,X_dt[,alls[s,]]))
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
      }
    }
  }
  exclude_lb=exclude_lb[!is.na(exclude_lb)]
  exclude_ub=exclude_ub[!is.na(exclude_ub)]
  exclude_intervals=interval_union(Intervals(cbind(exclude_lb,exclude_ub)))
  if(lb_max>-Inf){
    exclude_intervals=interval_union(exclude_intervals,Intervals(c(-Inf,lb_max)))
  }
  if(ub_min<Inf){
    exclude_intervals=interval_union(exclude_intervals,Intervals(c(ub_min, Inf)))
  }
  
  se=sqrt(sum(eta^2))*sigmahat    
  value=c(t(eta)%*%y)
  
  f1=function(mu){return (cdfF(mu, se, df=n-sum(selected)-1, exclude_intervals,value)-alpha/2)}
  f2=function(mu){return (cdfF(mu, se, df=n-sum(selected)-1, exclude_intervals,value)-1+alpha/2)}
  L=nleqslv(value+qnorm(alpha/2)*se,f2)$x
  U=nleqslv(value+qnorm(1-alpha/2)*se, f1)$x
  return(c(min(L,U), max(L,U)))
}

# helper func: cdf F(t) of truncated normal, given mean=mu, sd=se, and domain intervals
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

# helper func: find the quadratic root
quad <- function(A,B,C)
{
  answer <- c((-B - sqrt(B^2 - 4 * A * C)) / (2 * A),
              (-B + sqrt(B^2 - 4 * A * C)) / (2 * A))
  return(answer)
}

# Function: generate the 1st-order auto-regressive correlation matrix
ar1_matrix <- function(p, rho) {
  #input:
  #- p: number of variables
  #- rho: correlation paramter
  
  #output:
  #- a 1st-order auto-regressive correlation matrix
  
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p-1))
  rho^exponent
}

# Function generates an equal correlation matrix.
equal_cor_matrix <- function(p, rho) {
  #input:
  #- p: number of variables
  #- rho: correlation paramter

  #output:
  #- an equi-correlation matrix
  
  rho*matrix(rep(1, p ^ 2), nrow = p, ncol = p) + (1-rho)*diag(p)
}

# Function: Get standard error of mean response for new data samples, using ground truth sigma=1
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
