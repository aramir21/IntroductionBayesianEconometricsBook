remove(list = ls())
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/MateoV1/BS-UId/Data")
set.seed(12345)
n<-1000 
simnegbin =
  function(X, beta, alpha) {
    # Simulate from the Negative Binomial Regression
    lambda = exp(X %*% beta)
    y=NULL
    for (j in 1:length(lambda))
      y = c(y,rnbinom(1,mu = lambda[j],size = alpha))
    return(y)
  }
nvar=3 
alpha = 5
Vbeta = diag(nvar)*0.01
# Construct the regdata (containing X)
simnegbindata = NULL
beta = c(0.6,0.2,-0.5)
X = cbind(rep(1,n),rnorm(n,mean=2,sd=0.5),rnorm(n,mean=1,sd=1.5))
#simnegbindata = list(y=simnegbin(X,beta,alpha), X=X, beta=beta)
y=simnegbin(X,beta,alpha)
SimNegBin<-cbind(y,X[,-1])
colnames(SimNegBin)<-c("y","x1","x2")
write.csv(SimNegBin,file="17SimNegBinmodel.csv", row.names=FALSE)
library(MASS)
summary(m1 <- glm.nb(y ~ x1 + x2, data = as.data.frame(SimNegBin)))

