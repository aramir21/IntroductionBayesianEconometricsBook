remove(list = ls())
setwd("C:/Users/aramir21/Desktop/GUI/DataSim")
set.seed(12345)
n=100
X=cbind(rep(1,n),matrix(runif(3*n),n,3)); beta=c(1,2,-1,1.5); sigsq=.25
y=X%*%beta+rnorm(n,sd=sqrt(sigsq))
no<-26
Xother<-matrix(rnorm(n*no),n,no)
SimNormal<-as.data.frame(cbind(y,X[,-1],Xother))
SimNormal[1:10,]
write.csv(SimNormal,file="512SimNormalBMA.csv", row.names=FALSE)


