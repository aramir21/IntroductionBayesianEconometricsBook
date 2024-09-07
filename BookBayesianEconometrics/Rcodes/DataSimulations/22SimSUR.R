setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/BS-UId/Multivariate")
library(bayesm)
rm(list=ls())
set.seed(66)
beta1=c(1,2)
beta2=c(1,-1,-2)
nobs=1000
nreg=2
iota=c(rep(1,nobs))
X1=cbind(iota,runif(nobs))
X2=cbind(iota,runif(nobs),runif(nobs))
Sigma=matrix(c(.5,.2,.2,.5),ncol=2)
U=chol(Sigma)
E=matrix(rnorm(2*nobs),ncol=2)%*%U
y1=X1%*%beta1+E[,1]
y2=X2%*%beta2+E[,2]
dat<-cbind(y1,y2,X1,X2)
write.csv(dat,file="22SimSUR.csv")
