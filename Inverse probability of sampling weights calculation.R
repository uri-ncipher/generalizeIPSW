#Function that computes the IPSW estimator, its variance and 95% CI 
#data: combined trial (already subsetted so outcome data complete) 
#and cohort data with indicator for trial participation (s), covariates, outcome and predictors 
#nt: size of target population 
#selvar: vector of variables for sampling score model in combined trial and cohort data 
#trt: treatment or exposure variable in the trial 
#outcome: outcome variable in the trial 
ipsw.fun <- function(data,nt,selvar,trt,outcome) { 
  cohort=data[which(data$s==0),] 
  m=dim(cohort)[1] 
  trial=data[which(data$s==1),] 
  n=dim(trial)[1] 
  #define weights for selection propensity score model 
  trial$pw<-1 
  cohort$pw<-m/(nt-n) 
  #combine trial and cohort (S,Z) to estimate propensity scores 
  both<-rbind(cohort,trial) 
  b<-dim(both)[1] 
  #estimate selection propensity scores using logistic regression 
  #using data from cohort and trial stacked 
  mylogit <- glm(s ~ selvar, data = both, family = "binomial",weights=(both$pw)^(-1)) 
  both$p <- predict(mylogit,type = "response") 
  trial$p<-both$p[which(both$s==1)] 
  #Compute the IPSW estimator when S=1 
  #IPSW estimator: denominator 
  w1<-sum(trt/trial$p) 
  w2<-sum((1-trt)/trial$p) 
  #get mu(1) 
  trial$mu1i<- (trt*outcome)/trial$p 
  mu1hat <- (w1)^(-1)*sum(trial$mu1i) 
  #get mu(0) 
  trial$mu0i <- ((1-trt)*outcome)/trial$p 
  mu0hat<- (w2)^(-1)*sum(trial$mu0i) 
  PATEhat=mu1hat-mu0hat 
  #Estimated Variance for PATE when weights estimated  
  #Use sandwich estimator (Derived in Appendix A) 
  #estimate of derivative of weights 
  #first derivatives 
  #wrt beta1,beta10 
  #vector of covariates 
  both$int=rep(1,b) 
  trial$int=rep(1,n) 
  covb=as.matrix(cbind(both$int,selvar)) 
  covt<-as.matrix(cbind(trial$int,selvar[which(both$s==1),])) 
  both$wbeta <- covb*both$p*(1-both$p) 
  trial$wbeta<- covt*trial$p*(1-trial$p) 
  size=dim(covb)[2] 
  size2=size+2 
  bwbeta2<-list() 
  twbeta2<-list() 
  for(i in 1:size){ 
    bwbeta2[[i]]<-covb*(both$wbeta[,i]-2*both$wbeta[,i]*both$p) 
    twbeta2[[i]]<-covt*(trial$wbeta[,i]-2*trial$wbeta[,i]*trial$p) 
  } 
  Aw <- matrix(rep(NA,size*size),size,size) 
  for(i in 1:size){ 
    for(k in 1:size){ 
      Aw[i,k] <- sum(both$pw^(-1)*(-both$wbeta[,i]*both$wbeta[,k]*(both$p*(1-both$p) 
                                                                   + (both$s-both$p)*(1-2*both$p))/(both$p^2*(1-both$p)^2))) 
      + sum(both$pw^(-1)*(bwbeta2[[k]][,i]*(both$s-both$p))/(both$p*(1-both$p))) 
    } 
  }wAhat < -matrix(rep(NA,size2*size2),size2,size2) 
  wAhat[1,1]<-(b)^(-1)*sum((trt/trial$p)) 
  wAhat[1,2]<-0 
  wAhat[2,1]<-0 
  wAhat[2,2]<-(b)^(-1)*sum(((1-trt))/(trial$p)) 
  for(k in 1:size){ 
    wAhat[1,k+2] <-(b)^(-1)*sum((trt*(outcome-mu1hat)*trial$wbeta[,k])/trial$p^2) 
  }for(k in 1:size){ 
    wAhat[2,k+2] <-(b)^(-1)*sum(((1-trt)*(outcome-mu0hat)*trial$wbeta[,k])/trial$p^2) 
  }for(k in 3:size2){ 
    wAhat[k,1] <-0 
    wAhat[k,2] <-0 
  }for(i in 3:size2){ 
    for(k in 3:size2){ 
      wAhat[i,k]<-(b)^(-1)*(-Aw[i-2,k-2]) 
    } 
  }wBhat <- matrix(rep(NA,size2*size2),size2,size2) 
  wBhat[1,1]<-(b)^(-1)*sum((trt*(outcome-mu1hat)^2)/(trial$p)^2) 
  wBhat[1,2]<-wBhat[2,1]<-(b)^(-1)*sum((trt*(outcome-mu1hat))/trial$p*sum(((1-trt)*(outcome-mu0hat))/trial$p)) 
  wBhat[2,2]<-(b)^(-1)*sum(((1-trt)*(outcome-mu0hat)^2)/(trial$p)^2) 
  for(k in 1:size){ 
    wBhat[1,k+2]<-(b)^(-1)*sum((trt*(outcome-mu1hat)*trial$wbeta[,k])/(trial$p)^2) 
    wBhat[2,k+2]<-(b)^(-1)*sum(((1-trt)*(outcome-mu0hat)*trial$wbeta[,k])/(trial$p)^2) 
    wBhat[k+2,1]<-wBhat[1,k+2] 
    wBhat[k+2,2]<-wBhat[2,k+2] 
  }for(i in 1:size){ 
    for(k in 1:size){ 
      wBhat[i+2,k+2]<-(b)^(-1)*sum((((both$s-both$p)^2)/(both$p^2*(1-both$p)^2)) 
                                   *both$wbeta[,i]*both$wbeta[,k]*both$pw^(-2)) 
    } 
  }wsigmaM<-solve(wAhat)%*%wBhat%*%t(solve(wAhat)) 
  wsigma<-wsigmaM[1,1]+wsigmaM[2,2]-2*wsigmaM[2,1] 
  wse<-sqrt(wsigma)/sqrt(b) 
  PATE_LCL <-PATEhat-1.96*wse 
  PATE_UCL <-PATEhat+1.96*wse 
  result <-c(PATEhat, mu1hat,mu0hat,wsigma,wse, PATE_LCL, PATE_UCL) 
  names(result)<-c("PATEhat", "mu1hat","mu0hat","wsigma","wse", "PATE_LCL", "PATE_UCL") 
  return(result) 
} 