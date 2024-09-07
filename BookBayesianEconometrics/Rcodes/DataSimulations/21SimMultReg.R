setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/BS-UId/Multivariate")
library(bayesm)

set.seed(66)
n=1000
m=2
X=cbind(rep(1,n),runif(n),runif(n))
k=ncol(X)
B=matrix(c(1,2,-1,3,2,1),ncol=m)
Sigma=matrix(c(1,.5,.5,1),ncol=m); RSigma=chol(Sigma)
Y=X%*%B+matrix(rnorm(m*n),ncol=m)%*%RSigma
dat<-cbind(Y,X)
colnames(dat)<-c("y1","y2","x1","x2","x3")
write.csv(dat,file="21SimMultivariate.csv")
