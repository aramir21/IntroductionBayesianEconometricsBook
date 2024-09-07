#========================================
# Hierarchical Gaussian Linear Regression
#========================================
#== Generating data
# Constants
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/BEsmarterV11/Data")

library(MCMCpack)
nobs <- 1000
nspecies <- 20
species <- c(1:nspecies,sample(c(1:nspecies),(nobs-nspecies),replace=TRUE))
# Covariates
X1 <- runif(n=nobs,min=0,max=10)
X2 <- runif(n=nobs,min=0,max=10)
X <- cbind(rep(1,nobs),X1,X2)
W <- X
# Target parameters
# beta
beta.target <- matrix(c(0.1,0.3,0.2),ncol=1)
# Vb
Vb.target <- c(0.5,0.2,0.1)
# b
b.target <- cbind(rnorm(nspecies,mean=0,sd=sqrt(Vb.target[1])),
                  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[2])),
                  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[3])))
# sigma2
sigma2.target <- 0.02
# Response
Y <- vector()
for (n in 1:nobs) {
  Y[n] <- rnorm(n=1,
                mean=X[n,]%*%beta.target+W[n,]%*%b.target[species[n],],
                sd=sqrt(sigma2.target))
}
# Data-set
Data <- as.data.frame(cbind(Y,X1,X2,species))
plot(Data$X1,Data$Y)
#== Call to MCMChregress
model <- MCMChregress(fixed=Y~X1+X2, random=~X1+X2, group="species",
                      data=Data, burnin=1000, mcmc=1000, thin=1,verbose=1,
                      seed=NA, beta.start=0, sigma2.start=1,
                      Vb.start=1, mubeta=0, Vbeta=1.0E6,
                      r=3, R=diag(c(1,0.1,0.1)), nu=0.001, delta=0.001)
#== MCMC analysis
# Graphics
pdf("Posteriors-MCMChregress.pdf")
plot(model$mcmc)
dev.off()
# Summary
summary(model$mcmc)
# Predictive posterior mean for each observation
model$Y.pred
# Predicted-Observed
plot(Data$Y,model$Y.pred)
abline(a=0,b=1)
write.csv(Data,file="SimLogitudinalNormal.csv", row.names=FALSE)
## End(Not run)
