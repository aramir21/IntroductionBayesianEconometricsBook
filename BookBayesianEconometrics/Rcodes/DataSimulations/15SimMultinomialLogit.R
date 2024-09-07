#Maximum Likelihood Conditional Logit Model
remove(list = ls())
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/MateoV1/BS-UId/Data")
set.seed(12345)
MD<-1000  #Dataset Size
B<-c(0.5,0.8,-3)
B1<-c(-2.5,-3.5,0)
B2<-c(1,1,0)
set.seed(12345)
X1<-matrix(cbind(rnorm(MD,0,1),rnorm(MD,0,1),rnorm(MD,0,1)),MD,length(B))
X2<-matrix(cbind(rnorm(MD,0,1),rnorm(MD,0,1),rnorm(MD,0,1)),MD,length(B))
X3<-matrix(cbind(rnorm(MD,0,1),rnorm(MD,0,1),rnorm(MD,0,1)),MD,length(B))
X4<-matrix(rnorm(MD,1,1),MD,1)
V1<-B2[1]+X1%*%B+B1[1]*X4
V2<-B2[2]+X2%*%B+B1[2]*X4
V3<-B2[3]+X3%*%B+B1[3]*X4
suma<-exp(V1)+exp(V2)+exp(V3)
p1<-exp(V1)/suma
p2<-exp(V2)/suma
p3<-exp(V3)/suma
p<-cbind(p1,p2,p3)
y<- apply(p,1, function(x)sample(1:3, 1, prob = x, replace = TRUE))
table(y)
SimMultLogit<-cbind(y,X1[,1],X2[,1],X3[,1],X1[,2],X2[,2],X3[,2],X1[,3],X2[,3],X3[,3],X4)
colnames(SimMultLogit)<-c("y","x11","x12","x13","x21","x22","x23","x31","x32","x33","x4")
write.csv(SimMultLogit,file="15SimMultLogitmodel.csv", row.names=FALSE)
library(mlogit)
mode<-as.factor(y)
dat<-data.frame(mode,X1[,1],X2[,1],X3[,1],X1[,2],X2[,2],X3[,2],X1[,3],X2[,3],X3[,3],X4)
colnames(dat)<-c("mode","V1.1","V1.2","V1.3","V2.1","V2.2","V2.3","V3.1","V3.2","V3.3","V4")
attach(dat)
Exper<- mlogit.data(dat, shape = "wide", varying=2:10, choice = "mode")
res.mlogit<- mlogit(mode ~ V1 + V2 + V3 | V4, data=Exper, reflevel="3", probit=FALSE)
summary(res.mlogit)

Data<- read.table(file="15SimMultLogitmodel.csv",header=TRUE,sep=",")
Xa<-as.matrix(Data[,2:10])
Xd<-as.matrix(Data[,11])
XMPP<- createX(3, na=3, nd=1, Xa=Xa, Xd=Xd, INT = TRUE, DIFF = FALSE, base = 3)
#simout=simmnp(X,p,500,beta,Sigma)
Data1=list(y=as.vector(Data[,1]),X=XMPP,p=3)
Mcmc1=list(R=1000,keep=1)
out=rmnlIndepMetrop(Data=Data1,Mcmc=Mcmc1)
summary(out$betadraw)

