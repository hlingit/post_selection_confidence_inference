library(intervals)
library(nleqslv)

postICci=function(X_dt,y,selected, alls, new_xpoint, criteria='aic',alpha=0.05,sigmaKnown=F,sigmahat=NA){
  #  Input Parameters:
  #  selected: the AIC selected model (vector: boolean element indicating whether the predictor is included or not)
  #  alls: all other models to be compard with AIC-selected model (matrix: each row represent a model)
  #  sigmaKnown: if the noise level sigma^2 is known (boolean)
  #  sigmahat: the ground truth noise level if sigmaKnown=T; default NA
  
  #  Output:
  #  (1-alpha)-level confidence interval: c(L, U)
  
  n=length(y)
  if(criteria=='aic'){
    cons=2
  }else if(criteria=='bic'){
    cons=log(n)
  }else if(criteria=='aicc'){
    cons=2 #note: aicc=aic+(2k^2+2k)/(n-k-1)
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
    wt=exp(cons*(s_hat-s_tp)/n)
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
  
  f1=function(mu){return (cdfF(mu, se, df=n-sum(selected)-1, exclude_intervals,value, sigmaKnown = sigmaKnown)-alpha/2)}
  f2=function(mu){return (cdfF(mu, se, df=n-sum(selected)-1, exclude_intervals,value,sigmaKnown = sigmaKnown)-1+alpha/2)}
  L=nleqslv(value+qnorm(alpha/2)*se,f2)$x
  U=nleqslv(value+qnorm(1-alpha/2)*se, f1)$x
  return(c(min(L,U), max(L,U)))
}

#helper func
cdfF=function(mu, se, df=NA, exclude_intervals,value,sigmaKnown=F){
  if(!sigmaKnown){
    #rescale the truncated intervals so we can use truncated t
    rescale_excludeInt=(exclude_intervals@.Data-mu)/se
    rescale_value=(value-mu)/se
    normalizer=1-sum(pt(rescale_excludeInt[,2], df=df)-pt(rescale_excludeInt[,1],df=df))
    int_left_endpoints=pt(rescale_excludeInt[,1], df=df)
    if(length(int_left_endpoints)>1){
      int_left_endpoints=int_left_endpoints-c(0, pt(rescale_excludeInt[-length(int_left_endpoints),2], df=df))
      int_left_endpoints=cumsum(int_left_endpoints)
    }
    int_left_endpoints=int_left_endpoints/normalizer
    value_loc=which(rescale_excludeInt[,1]>rescale_value)[1]-1
    if(is.na(value_loc)){
      value_loc=length(int_left_endpoints)
    }
    result=int_left_endpoints[value_loc]+(pt(rescale_value,df=df)-pt(rescale_excludeInt[value_loc,2], df=df))/normalizer
  }else{
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
  }
  return(result)
}

#helper func: find the quadratic root
quad <- function(A,B,C)
{
  answer <- c((-B - sqrt(B^2 - 4 * A * C)) / (2 * A),
              (-B + sqrt(B^2 - 4 * A * C)) / (2 * A))
  return(answer)
}