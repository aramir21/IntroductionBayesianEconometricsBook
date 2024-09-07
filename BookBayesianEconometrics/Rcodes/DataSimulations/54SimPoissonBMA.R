remove(list = ls())
setwd("C:/Users/aramir21/Desktop/GUI/DataSim")
set.seed(12345)
n<-1000 
B<-c(0.5,0.8,-1.2)
X<-matrix(cbind(rep(1,n),rnorm(n,0,1),rnorm(n,0,1)),n,length(B))
y1<-exp(X%*%B)
y<-rpois(n,y1)
hist(y)
nXgar<-25
Xgar<-matrix(rnorm(nXgar*n),n,nXgar)
SimLogit<-as.data.frame(cbind(y,X[,-1],Xgar))
SimLogit[1:10,]
write.csv(SimLogit,file="54SimPoissonBMA.csv", row.names=FALSE)
reg<-glm(y~X[,-1]+Xgar,data=SimLogit,family=poisson(link = "log"))
summary(reg)
