#######Multivariate probit model#######
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/BS-UId/Multivariate")
remove(list = ls())
n<-5000  #Individuals
p<-2  #Number dependent variables (For instance to buy or not to buy different products)
xo<- 2  #Number of choice dependent variables (For instance price of products)
xi<- 1  #Number of regressors that dependt on individuals (For instance income)
B1<-c(0.5,-1.2,0.7,0.8)
B2<-c(1.5,-0.8,0.5,1.2)
Cor<-matrix(c(1,0.5,0.5,1),p,p)

yX<-NULL

for (i in 1:n){
  ei<-mvrnorm(1,c(0,0),Cor)
  xs<-rnorm(xi)
  x1i<-c(1,rnorm(xo),xs)
  yl1i<-sum(x1i*B1)
  y1i<-yl1i+ei[1]>0
  x2i<-c(1,rnorm(xo),xs)
  yl2i<-sum(x2i*B2)
  y2i<-yl2i+ei[2]>0
  yXi<-rbind(c(y1i,x1i),c(y2i,x2i))
  yXi<-cbind(i,yXi)
  colnames(yXi)<-c("id","y","cte","x1o","x2o","xi")
  yX<-rbind(yX,yXi)
}

write.csv(yX,file="24SimMultProbit.csv")
