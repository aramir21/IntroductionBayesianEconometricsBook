remove(list = ls())
setwd("C:/Users/aramir21/Desktop/GUI/DataSim")
set.seed(12345)
n<-1000 
B<-c(0.5,0.5,0.5)
X<-matrix(cbind(rep(1,n),rnorm(n,0,0.5),rnorm(n,0,0.5)),n,length(B))
y1<-exp(X%*%B)
y<-rgamma(n,y1,rate=1)
summary(y)
plot(density(y))
nXgar<-25
Xgar<-matrix(rnorm(nXgar*n),n,nXgar)
SimLogit<-as.data.frame(cbind(y,X[,-1],Xgar))
SimLogit[1:10,]
write.csv(SimLogit,file="53SimGammaBMA.csv", row.names=FALSE)
reg<-glm(y~X[,-1]+Xgar,data=SimLogit,family=Gamma(link = "inverse"))
summary(reg)
