remove(list = ls())
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/MateoV1/BS-UId/Data")
set.seed(12345)
n<-1000 
B<-c(0.5,0.8,-1.2)
X<-matrix(cbind(rep(1,n),rnorm(n,0,1),rnorm(n,0,1)),n,length(B))
u<-rlogis(n,0,1)
y1<-X%*%B+u
y<-y1>0
table(y)
SimLogit<-as.data.frame(cbind(y,X[,-1]))
names(SimLogit)<-c("y","x1","x2")
SimLogit[1:10,]
write.csv(SimLogit,file="12SimLogitmodel.csv", row.names=FALSE)
reg<-glm(y~x1+x2,data=SimLogit,family=binomial(link="logit"))
summary(reg)
