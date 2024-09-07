rm(list = ls())
set.seed(010101)
########################## Linear regression: Value of soccer players (heteroskedasticity) ##########################
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
v <- 5
# MCMC parameters
mcmc <- 5000
burnin <- 5000
tot <- mcmc + burnin
thin <- 1
# Posterior distributions programming the Gibbs sampling
# Auxiliary parameters
an <- a0 + N
v1n <- v + 1
# Gibbs sampling functions
PostSig2 <- function(Beta, tau){
  dn <- d0 + t(y - X%*%Beta)%*%diag(tau)%*%(y - X%*%Beta)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBeta <- function(sig2, tau){
  Bn <- solve(B0i + sig2^(-1)*t(X)%*%diag(tau)%*%X)
  bn <- Bn%*%(B0i%*%b0 + sig2^(-1)*t(X)%*%diag(tau)%*%y)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostTau <- function(sig2, Beta, i){
  v2n <- v + sig2^(-1)*(y[i]-X[i,]%*%Beta)^2
  taui <- rgamma(1, v1n, v2n)
  return(taui)
}
PostBetas <- matrix(0, mcmc+burnin, k)
PostSigma2 <- rep(0, mcmc+burnin)
Beta <- rep(0, k)
tau <- rep(1, N)
for(s in 1:tot){
  sig2 <- PostSig2(Beta = Beta, tau = tau)
  PostSigma2[s] <- sig2
  Beta <- PostBeta(sig2 = sig2, tau = tau)
  PostBetas[s,] <- Beta
  tau <- sapply(1:N, function(i){PostTau(sig2 = sig2, Beta = Beta, i)})
}
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
colnames(PosteriorBetas) <- c("Intercept", "Perf", "Age", "Age2", "NatTeam", "Goals", "Exp", "Exp2")
summary(coda::mcmc(PosteriorBetas))
PosteriorSigma2 <- PostSigma2[keep]
summary(coda::mcmc(PosteriorSigma2))
summary(coda::mcmc(exp(PosteriorBetas[,5])-1))

########################## Probit model: Determinants of hospitalization ##########################
rm(list = ls())
set.seed(010101)
Data <- read.csv("DataApplications/2HealthMed.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(Data)
str(Data)
y <- Hosp # Dependent variables
X <- cbind(1, SHI, Female, Age, Age2, Est2, Est3, Fair, Good, Excellent) # Regressors
K <- dim(X)[2] # Number of regressors
N <- dim(X)[1]
# Hyperparameters
b0 <- rep(0, K) # Prio mean
B0 <- diag(K) # Prior precision (inverse of covariance)
B0i <- solve(B0)
mcmc <- 1000; burnin <- 500; thin <- 2; tot <- mcmc + burnin; keep <- seq(burnin, tot, thin)
# Posterior distributions programming the Gibbs sampling
# Auxiliary parameters
XtX <- t(X)%*%X
# Gibbs sampling functions
PostBeta <- function(Yl){
  Bn <- solve(B0i + XtX)
  bn <- Bn%*%(B0i%*%b0 + t(X)%*%Yl)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostYl <- function(Beta, i){
  Ylmean <- X[i,]%*%Beta
  if(y[i] == 1){
    Yli <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = Ylmean, sd = 1)
  }else{
    Yli <- truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = Ylmean, sd = 1)
  }
  return(Yli)
}
PostBetas <- matrix(0, mcmc+burnin, K)
Beta <- rep(0, K)
# create progress bar in case that you want to see iterations progress
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Yl <- sapply(1:N, function(i){PostYl(Beta = Beta, i)})
  Beta <- PostBeta(Yl = Yl)
  PostBetas[s,] <- Beta
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0), "% done"))

}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
colnames(PosteriorBetas) <- c("Intercept", "SHI", "Female", "Age", "Age2", "Est2", "Est3", "Fair", "Good", "Excellent")
summary(coda::mcmc(PosteriorBetas))


