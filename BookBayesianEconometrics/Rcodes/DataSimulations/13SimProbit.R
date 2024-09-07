remove(list = ls())
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/MateoV1/BS-UId/Data")
set.seed(12345)
n<-1000 
B<-c(0.5,0.8,-1.2)
X<-matrix(cbind(rep(1,n),rnorm(n,0,1),rnorm(n,0,1)),n,length(B))
u<-rnorm(n)
y1<-X%*%B+u
y<-y1>0
table(y)
SimProbit<-as.data.frame(cbind(y,X[,-1]))
names(SimProbit)<-c("y","x1","x2")
SimProbit[1:10,]
write.csv(SimProbit,file="12SimProbitmodel.csv", row.names=FALSE)
glm(y~x1+x2,data=SimProbit,family=binomial(link="probit"))
