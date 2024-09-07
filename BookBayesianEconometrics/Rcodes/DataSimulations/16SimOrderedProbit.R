remove(list = ls())
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/MateoV1/BS-UId/Data")
set.seed(12345)
n<-1000 
B<-c(0.5,0.8,-1.2)
X<-matrix(cbind(rnorm(n,0,1),rnorm(n,0,1),rnorm(n,0,1)),n,length(B))
u<-rnorm(n,0,1)
Y<-X%*%B+u
a1<- 0
a2<- 1.3
y<-ifelse(Y<=a1,1,ifelse(Y>a2,3,2))
table(y)
SimOrdered<-as.data.frame(cbind(y,X))
names(SimOrdered)<-c("y","x1","x2","x3")
SimOrdered[1:10,]
write.csv(SimOrdered,file="16SimOrderedmodel.csv", row.names=FALSE)
polr(as.factor(y)~X, method = "probit")
