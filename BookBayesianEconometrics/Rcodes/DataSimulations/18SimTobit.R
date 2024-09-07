remove(list = ls())
library(AER)
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/MateoV1/BS-UId/Data")
set.seed(12345)
n=1000
X=cbind(rep(1,n),matrix(runif(3*n),n,3)); beta=c(1,2,-1,1.5); sigsq=.25
y=X%*%beta+rnorm(n,sd=sqrt(sigsq))
summary(y)
y.ast<-replace(y, y<0, 0)
y.ast<-replace(y.ast, y.ast>3, 3)
y[y<0]<- 0
y[y>3]<- 3
cbind(y.ast,y)
SimTobit<-as.data.frame(cbind(y.ast,X[,-1]))
names(SimTobit)<-c("y","x1","x2","x3")
SimTobit[1:10,]
write.csv(SimTobit,file="18SimTobitmodel.csv", row.names=FALSE)
tobit<-tobit(y~x1+x2+x3,data=SimTobit,left=0, right=Inf)
summary(tobit)

