## Not run:
#========================================
# Hierarchical Poisson Linear Regression
#========================================
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/BEsmarterV11/Data")
library(MCMCpack)

#== Generating data
# Constants
nobs <- 1000
nspecies <- 20
species <- c(1:nspecies,sample(c(1:nspecies),(nobs-nspecies),replace=TRUE))
# Covariates
X1 <- runif(n=nobs,min=-1,max=1)
X2 <- runif(n=nobs,min=-1,max=1)
X <- cbind(rep(1,nobs),X1,X2)
W <- X
# Target parameters
# beta
beta.target <- matrix(c(0.1,0.1,0.1),ncol=1)
# Vb
Vb.target <- c(0.05,0.05,0.05)
# b
b.target <- cbind(rnorm(nspecies,mean=0,sd=sqrt(Vb.target[1])),
                  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[2])),
                  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[3])))
# Response
lambda <- vector()
Y <- vector()
for (n in 1:nobs) {
  lambda[n] <- exp(X[n,]%*%beta.target+W[n,]%*%b.target[species[n],])
  Y[n] <- rpois(1,lambda[n])
}
# Data-set
Data <- as.data.frame(cbind(Y,lambda,X1,X2,species))
plot(Data$X1,Data$lambda)
#== Call to MCMChpoisson
model <- MCMChpoisson(fixed=Y~X1+X2, random=~X1+X2, group="species",
                      data=Data, burnin=5000, mcmc=1000, thin=1,verbose=1,
                      seed=NA, beta.start=0, sigma2.start=1,
                      Vb.start=1, mubeta=0, Vbeta=1.0E6,
                      r=3, R=diag(c(0.1,0.1,0.1)), nu=0.001, delta=0.001, FixOD=0)
#== MCMC analysis
# Graphics
pdf("Posteriors-MCMChpoisson.pdf")
plot(model$mcmc)
dev.off()
# Summary
summary(model$mcmc)
# Predictive posterior mean for each observation
model$lambda.pred
# Predicted-Observed
plot(Data$lambda,model$lambda.pred)
abline(a=0,b=1)
Data <- as.data.frame(cbind(Y,X1,X2,species))
write.csv(Data,file="SimLogitudinalPoisson.csv", row.names=FALSE)
## #Not run
## #You can also compare with lme4 results
## #== lme4 resolution
## library(lme4)
## model.lme4 <- lmer(Y~X1+X2+(1+X1+X2|species),data=Data,family="poisson")
## summary(model.lme4)
## plot(fitted(model.lme4),model$lambda.pred,main="MCMChpoisson/lme4")
## abline(a=0,b=1)
## End(Not run)
