rm(list = ls())
set.seed(010101)
########################## Linear regression: Value of soccer players ##########################
Data <- read.csv("DataApplications/1ValueFootballPlayers.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(Data)
y <- log(Value) 
# Value: Market value in Euros (2017) of soccer players
# Regressors quantity including intercept
X <- cbind(1, Perf, Age, Age2, NatTeam, Goals, Exp, Exp2)
# Perf: Performance. Perf2: Performance squared. Age: Age; Age: Age squared. 
# NatTeam: Indicator of national team. Goals: Scored goals. Goals2: Scored goals squared
# Exp: Years of experience. Exp2: Years of experience squared. Assists: Number of assists
k <- dim(X)[2]
N <- dim(X)[1]
# Hyperparameters
d0 <- 0.001/2
a0 <- 0.001/2
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
# MCMC parameters
mcmc <- 5000
burnin <- 5000
tot <- mcmc + burnin
thin <- 1
# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
posterior  <- MCMCpack::MCMCregress(y~X-1, b0=b0, B0 = B0i, c0 = a0, d0 = d0, burnin = burnin, mcmc = mcmc, thin = thin)
summary(coda::mcmc(posterior))
# Posterior distributions programming the Gibbs sampling
# Auxiliary parameters
XtX <- t(X)%*%X
bhat <- solve(XtX)%*%t(X)%*%y
an <- a0 + N
# Gibbs sampling functions
PostSig2 <- function(Beta){
  dn <- d0 + t(y - X%*%Beta)%*%(y - X%*%Beta)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBeta <- function(sig2){
  Bn <- solve(B0i + sig2^(-1)*XtX)
  bn <- Bn%*%(B0i%*%b0 + sig2^(-1)*XtX%*%bhat)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostBetas <- matrix(0, mcmc+burnin, k)
PostSigma2 <- rep(0, mcmc+burnin)
Beta <- rep(0, k)
for(s in 1:tot){
  sig2 <- PostSig2(Beta = Beta)
  PostSigma2[s] <- sig2
  Beta <- PostBeta(sig2 = sig2)
  PostBetas[s,] <- Beta
}
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
colnames(PosteriorBetas) <- c("Intercept", "Perf", "Age", "Age2", "NatTeam", "Goals", "Exp", "Exp2")
summary(coda::mcmc(PosteriorBetas))
PosteriorSigma2 <- PostSigma2[keep]
summary(coda::mcmc(PosteriorSigma2))
summary(coda::mcmc(exp(PosteriorBetas[,5])-1))

########################## Logit: Simulation ##########################
# Simulate data
rm(list = ls())
set.seed(010101)
N <- 10000 # Sample size
B <- c(0.5, 0.8, -1.2) # Population location parameters
x2 <- rnorm(N) # Regressor
x3 <- rnorm(N) # Regressor
X <- cbind(1, x2, x3) # Regressors
XB <- X%*%B
PY <- exp(XB)/(1 + exp(XB)) # Probability of Y = 1
Y <- rbinom(N, 1, PY) # Draw Y's
table(Y) # Frequency
# write.csv(cbind(Y, x2, x3), file = "DataSimulations/LogitSim.csv") # Export data

# MCMC parameters
iter <- 5000; burnin <- 1000; thin <- 5; tune <- 1
# Hyperparameters
K <- dim(X)[2] 
b0 <- rep(0, K)
c0 <- 1000
B0 <- c0*diag(K)
B0i <- solve(B0)
# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
RegLog <- MCMCpack::MCMClogit(Y~X-1, mcmc = iter, burnin = burnin, thin = thin, b0 = b0, B0 = B0i, tune = tune)
summary(RegLog)
# Posterior distributions programming the Metropolis-Hastings algorithm
# This function sets the M-H sampler using as default a hyperparameter mean equal to 0
# and a covariance equal to 1000 times a identity matrix, a tuning parameter equal to 1,
# 1000 post burn-in iterations, and the latter is equal to 500.
MHfunc <- function(y, X, b0 = rep(0, dim(X)[2] + 1), B0 = 1000*diag(dim(X)[2] + 1), tau = 1, 
                   iter = 6000, burnin = 1000, thin = 5){
  Xm <- cbind(1, X) # Regressors
  K <- dim(Xm)[2] # Number of location parameters
  BETAS <- matrix(0, iter + burnin, K) # Space for posterior chains
  Reg <- glm(y ~ Xm - 1, family = binomial(link = "logit")) # Maximum likelihood estimation
  BETA <- Reg$coefficients # Maximum likelihood parameter estimates 
  tot <- iter + burnin # Total iterations M-H algorithm
  COV <- vcov(Reg) # Maximum likelihood covariance matrix
  COVt <- tau^2*solve(solve(B0) + solve(COV)) # Covariance matrix for the proposal distribution
  Accep <- rep(0, tot) # Space for calculating the acceptance rate
  # create progress bar in case that you want to see iterations progress
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = tot, width = 300)
  for(it in 1:tot){
    BETAc <- BETA + MASS::mvrnorm(n = 1, mu = rep(0, K), Sigma = COVt) # Candidate location parameter
    likecand <- sum((Xm%*%BETAc) * Y - apply(Xm%*%BETAc, 1, function(x) log(1 + exp(x)))) # Log likelihood for the candidate
    likepast <- sum((Xm%*%BETA) * Y - apply((Xm%*%BETA), 1, function(x) log(1 + exp(x)))) # Log likelihood for the actual draw
    priorcand <- (-1/2)*crossprod((BETAc - b0), solve(B0))%*%(BETAc - b0) # Log prior for candidate
    priorpast <- (-1/2)*crossprod((BETA - b0), solve(B0))%*%(BETA - b0) # Log prior for actual draw
    alpha <- min(1, exp(likecand + priorcand - likepast - priorpast)) #Probability of selecting candidate
    u <- runif(1) # Decision rule for selecting candidate
    if(u < alpha){
      BETA <- BETAc # Changing reference for candidate if selected
      Accep[it] <- 1 # Indicator if the candidate is accepted
    } 
    BETAS[it, ] <- BETA # Saving draws
    setWinProgressBar(pb, it, title=paste( round(it/tot*100, 0),
                                          "% done"))
  }
  close(pb)
  keep <- seq(burnin, tot, thin)
  return(list(Bs = BETAS[keep[-1], ], AceptRate = mean(Accep[keep[-1]])))
}
Posterior <- MHfunc(y = Y, X = cbind(x2, x3), iter = iter, burnin = burnin, thin = thin) # Running our M-H function changing some default parameters.
paste("Acceptance rate equal to", round(Posterior$AceptRate, 2), sep = " ")
PostPar <- coda::mcmc(Posterior$Bs)
# Names
colnames(PostPar) <- c("Cte", "x1", "x2")
# Summary posterior draws
summary(PostPar)
# Trace and density plots
plot(PostPar)
# Autocorrelation plots
coda::autocorr.plot(PostPar)
# Convergence diagnostics
coda::geweke.diag(PostPar)
coda::raftery.diag(PostPar,q=0.5,r=0.05,s = 0.95)
coda::heidel.diag(PostPar)

########################## Probit: Hospitalizationn ##########################
rm(list = ls())
set.seed(010101)
mydata <- read.csv("DataApplications/2HealthMed.csv", header = T, sep = ",")
attach(mydata)
str(mydata)
y <- Hosp # Dependent variables
X <- cbind(1, SHI, Female, Age, Age2, Est2, Est3, Fair, Good, Excellent) # Regressors
K <- 10 # Number of regressors
# Hyperparameters
b0 <- rep(0, K) # Prio mean
B0i <- diag(K) # Prior precision (inverse of covariance)
Prior <- list(betabar = b0, A = B0i) # Prior list
Data <- list(y = y, X = X) # Data list
Mcmc <- list(R = 10000, keep = 1, nprint = 0) # MCMC parameters
RegProb <- bayesm::rbprobitGibbs(Data = Data, Prior = Prior, Mcmc = Mcmc) # Inference using bayesm package
PostPar <- coda::mcmc(RegProb$betadraw) # Posterior draws
colnames(PostPar) <- c("Cte", "SHI", "Female", "Age", "Age2", "Est2",
                       "Est3", "Fair", "Good", "Excellent") # Names
summary(PostPar) # Posterior summary

