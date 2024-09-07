## Not run:
#========================================
# Hierarchical Binomial Linear Regression
#========================================
setwd("C:/ANDRES/Portatil/Bayesian Econometrics/UserInterface/BEsmarterV11/Data")
library(MCMCpack)
#== inv.logit function
inv.logit <- function(x, min=0, max=1) {
  p <- exp(x)/(1+exp(x))
  p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
  return(p*(max-min)+min)
}
#== Generating data
# Constants
nobs <- 1000
nspecies <- 20
species <- c(1:nspecies,sample(c(1:nspecies),(nobs-nspecies),replace=TRUE))
# Covariates
X1 <- runif(n=nobs,min=-10,max=10)
X2 <- runif(n=nobs,min=-10,max=10)
X <- cbind(rep(1,nobs),X1,X2)
W <- X
# Target parameters
# beta
beta.target <- matrix(c(0.3,0.2,0.1),ncol=1)
# Vb
Vb.target <- c(0.5,0.05,0.05)
# b
b.target <- cbind(rnorm(nspecies,mean=0,sd=sqrt(Vb.target[1])),
                  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[2])),
                  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[3])))
# Response
theta <- vector()
Y <- vector()
for (n in 1:nobs) {
  theta[n] <- inv.logit(X[n,]%*%beta.target+W[n,]%*%b.target[species[n],])
  Y[n] <- rbinom(n=1,size=1,prob=theta[n])
}
# Data-set
Data <- as.data.frame(cbind(Y,theta,X1,X2,species))
plot(Data$X1,Data$theta)
#== Call to MCMChlogit
model <- MCMChlogit(fixed=Y~X1+X2, random=~X1+X2, group="species",
                    data=Data, burnin=5000, mcmc=1000, thin=1,verbose=1,
                    seed=NA, beta.start=0, sigma2.start=1,
                    Vb.start=1, mubeta=0, Vbeta=1.0E6,
                    r=3, R=diag(c(1,0.1,0.1)), nu=0.001, delta=0.001, FixOD=1)
#== MCMC analysis
# Graphics
pdf("Posteriors-MCMChlogit.pdf")
plot(model$mcmc)
dev.off()
# Summary
summary(model$mcmc)
# Predictive posterior mean for each observation
model$theta.pred
# Predicted-Observed
plot(Data$theta,model$theta.pred)
abline(a=0,b=1)
Data <- as.data.frame(cbind(Y,X1,X2,species))
write.csv(Data,file="SimLogitudinalLogit.csv", row.names=FALSE)
## #Not run
## #You can also compare with lme4 results
## #== lme4 resolution
## library(lme4)
## model.lme4 <- lmer(Y~X1+X2+(1+X1+X2|species),data=Data,family="binomial")
## summary(model.lme4)
## plot(fitted(model.lme4),model$theta.pred,main="MCMChlogit/lme4")
## abline(a=0,b=1)
## End(Not run)