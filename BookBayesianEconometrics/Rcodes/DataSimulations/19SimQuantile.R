remove(list = ls())
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/MateoV1/BS-UId/Data")
set.seed(12345)
n=1000
X=cbind(rep(1,n),matrix(runif(3*n),n,3)); beta=c(1,2,-1,1.5); sigsq=.25
y=X%*%beta+rnorm(n,sd=sqrt(sigsq))
SimNormal<-as.data.frame(cbind(y,X[,-1]))
names(SimNormal)<-c("y","x1","x2","x3")
SimNormal[1:10,]
write.csv(SimNormal,file="19SimQuantilemodel.csv", row.names=FALSE)
glm(y~x1+x2+x3,data=SimNormal,family="gaussian")
