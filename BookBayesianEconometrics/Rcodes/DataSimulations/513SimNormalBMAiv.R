setwd("C:/Users/aramir21/Desktop/GUI/DataSim")
set.seed(66)
simIV = function(delta1,delta2,beta0,betas1,betas2,beta2,Sigma,n,z) {
  eps = matrix(rnorm(3*n),ncol=3) %*% chol(Sigma)
  xs1 = z%*%delta1 + eps[,1]
  xs2 = z%*%delta2 + eps[,2]
  x2 = rnorm(dim(z)[1])
  y = beta0+betas1*xs1+betas2*xs2+beta2*x2 + eps[,3]
  X = as.matrix(cbind(xs1,xs2,1,x2)) 
  colnames(X)=c("x1en","x2en","cte","xex")
  y=matrix(y,dim(z)[1],1)
  colnames(y)=c("y")
  list(X=X,y=y)}
n = 1000 ; p=3 
z = matrix(runif(n*p),ncol=p)
rho31=.8; rho32=.5;
Sigma = matrix(c(1,0,rho31,0,1,rho32,rho31,rho32,1),ncol=3)
delta1 = c(4,-1,2); delta2=c(-2,3,-1); betas1 = .5; betas2 = -1; beta2 = 1; beta0=2
simiv = simIV(delta1,delta2,beta0,betas1,betas2,beta2,Sigma,n,z)
nW<- 18
W<- matrix(rnorm(nW*dim(z)[1]),dim(z)[1],nW)
YXW<-cbind(simiv$y,simiv$X,W)
write.csv(YXW,file='513SimNormalBMAivYXWNew.csv', row.names=FALSE)
Z<-z
write.csv(Z,file="513SimNormalBMAivZNew.csv", row.names=FALSE)


