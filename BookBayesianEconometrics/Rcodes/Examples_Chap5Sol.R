########################## Math test example ########################## 
rm(list = ls())
set.seed(010101)
N <- 50
# Sample size
muhat <- 102
# Sample mean
sig2hat <- 100
# Sample variance

# Hyperparameters
mu0 <- 100
sig20 <- 100
delta0 <- 0.001
alpha0 <- 0.001

MCMC <- 10000; burnin <- 1000; S <- MCMC + burnin
keep <- (burnin+1):S
# Posterior draws
alphan <- alpha0 + N
sig2Post <- rep(NA, S)
muPost <- rep(NA, S)
sig2 <- sig20
for(s in 1:S){
  sig2n <- (1/(sig2/N)+1/sig20)^(-1)
  mun <- sig2n*(muhat/(sig2/N)+mu0/sig20)
  mu <- rnorm(1, mun, sig2n^0.5)
  muPost[s] <- mu
  deltan <- N*(sig2hat + (muhat - mu)^2)
  sig2 <- invgamma::rinvgamma(1, shape = alphan, rate = deltan)
  sig2Post[s] <- sig2
}
sig2s <- coda::mcmc(sig2Post[keep]) 
mus <- coda::mcmc(muPost[keep]) 
summary(sig2s)
summary(mus)
hist(mus, main = "Histogram: Posterior mean", xlab = "Posterior mean", col = "blue", breaks = 50)
muPost_tq <- quantile(mus, c(0.025, 0.5, 0.975))
muPost_tq
cutoff <- 103
PmuPost_tcutoff <- mean(mus > cutoff)
PmuPost_tcutoff

########################## Metropolis-Hastings: Cauchy distribution ########################## 
rm(list = ls())
set.seed(010101)
S <- 100000
df <- 5
theta1 <- runif(S)
acept1 <- rep(0, S)
theta2 <- runif(S)
acept2 <- rep(0, S)
for (s in 2:S){
  thetac1 <- rnorm(1)
  thetac2 <- rt(1, df = df)
  a1 <- (dcauchy(thetac1)*dnorm(theta1[s-1]))/(dcauchy(theta1[s-1])*dnorm(thetac1))
  a2 <- (dcauchy(thetac2)*dt(theta2[s-1], df = df))/(dcauchy(theta2[s-1])*dt(thetac2, df = df))
  U1 <- runif(1)
  if(U1 <= a1){
    theta1[s] <- thetac1
    acept1[s] <- 1
  }else{
    theta1[s] <- theta1[s-1]
    acept1[s] <- 0
  }
  if(U1 <= a2){
    theta2[s] <- thetac2
    acept2[s] <- 1
  }else{

    theta2[s] <- theta2[s-1]
    acept2[s] <- 0
  }
}
mean(acept1); mean(acept2)
mean(theta1); sd(theta1)
mean(theta2); sd(theta2)
plot(coda::mcmc(theta1)); coda::autocorr.plot(coda::mcmc(theta1))
plot(coda::mcmc(theta2)); coda::autocorr.plot(coda::mcmc(theta2))
h <- hist(theta1, breaks=50, col="blue", xlab="x", main="Cauchy draws from a Metropolis-Hastings algorithm: Normal standard proposal")
pfit <- seq(min(theta1),max(theta1),length=50)
yfit<-dcauchy(pfit)
yfit <- yfit*diff(h$mids[1:2])*length(theta1)
lines(pfit, yfit, col="red", lwd=2)

h <- hist(theta2, breaks=50, col="blue", xlab="x", main="Cauchy draws from a Metropolis-Hastings algorithm: Student's t proposal")
pfit <- seq(min(theta2),max(theta2),length=50)
yfit<-dcauchy(pfit)
yfit <- yfit*diff(h$mids[1:2])*length(theta2)
lines(pfit, yfit, col="red", lwd=2)


########################## Hamiltonian Monte Carlo: Bimodal distribution ########################## 
# https://hedibert.org/wp-content/uploads/2021/02/hmc-example.html
rm(list = ls()); set.seed(010101)
# Posterior distribution
PostDist <- function(theta){
  theta1 <- theta[1]; theta2 <- theta[2]
  post <- exp(-0.5*(theta1^2*theta2^2+theta1^2+theta2^2-8*theta1-8*theta2))
  return(post)
}
# Negative of the posterior to calculate the covariance matrix
LogPostDistNeg <- function(theta){
  theta1 <- theta[1]; theta2 <- theta[2]
  post <- -0.5*(theta1^2*theta2^2+theta1^2+theta2^2-8*theta1-8*theta2)
  return(-post)
}
theta0 <- c(0,0)
ResOpt <- optim(theta0, LogPostDistNeg, hessian = TRUE)
VarPost <- solve(ResOpt$hessian)
# Metropolis-Hastings
MH <- function(theta, c2){
  thetac <- MASS::mvrnorm(1, mu = theta, Sigma = c2*diag(2))
  a <- PostDist(thetac)/PostDist(theta)
  U <- runif(1)
  if(U <= a){
    theta <- thetac
    accept <- 1
  }else{
    theta <- theta
    accept <- 0
  }
  return(list(theta = theta, accept = accept))
}
# Posterior draws M-H
S <- 20000; burnin <- 1000; thin <- 20; tot <- S + burnin
K <- 2; c2 <- 1.5
thetaPostMH <- matrix(NA, tot, K)
AcceptMH <- rep(NA, tot)
thetaMH <- c(0, 0)
for(s in 1:tot){
  ResMH <- MH(thetaMH, c2 = c2)
  thetaMH <- ResMH$theta
  thetaPostMH[s,] <- thetaMH
  AcceptMH[s] <- ResMH$accept
}
keep <- seq((burnin), tot, thin)
mean(AcceptMH[keep])
thetaPostMHMCMC <- coda::mcmc(thetaPostMH[keep,])
plot(thetaPostMHMCMC)
coda::autocorr.plot(thetaPostMHMCMC)
# Contour plot
ngrid <- 400
theta1 <- seq(-1, 6, length = ngrid)
theta2 <- seq(-1, 6, length = ngrid)
f <- matrix(0, ngrid, ngrid)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    f[i,j] = PostDist(c(theta1[i],theta2[j]))
  }
}
plot(thetaPostMH[keep,], xlim=range(theta1), ylim=range(theta2), pch=16, col=grey(0.8), xlab=expression(theta[1]), ylab=expression(theta[2]))
contour(theta1,theta2,f, drawlabels=FALSE, add=TRUE, col = "blue", lwd = 1.2)
title("Random walk Metropolis-Hastings")
# Hamiltonian Monte Carlo
HMC <- function(theta, epsilon, L, M){
  # L <- ceiling(1/epsilon)
  Minv <- solve(M); thetat <- theta
  K <- length(thetat)
  mom <- t(mvtnorm::rmvnorm(1, rep(0, K), M))
  logPost_Mom_t <- log(PostDist(thetat)) +  mvtnorm::dmvnorm(t(mom), rep(0, K), M, log = TRUE)  
  theta1 <- theta[1]; theta2 <- theta[2] 
  for(l in 1:L){
    if(l == 1 | l == L){
      mom <- mom + 0.5*epsilon*(-0.5*c(2*theta1*theta2^2+2*theta1-8, 2*theta2*theta1^2+2*theta2-8))
      theta <- theta + epsilon*Minv%*%mom
    }else{
      mom <- mom + epsilon*(-0.5*c(2*theta1*theta2^2+2*theta1-8, 2*theta2*theta1^2+2*theta2-8))
      theta <- theta + epsilon*Minv%*%mom
    }
  }
  logPost_Mom_star <- log(PostDist(theta)) +  mvtnorm::dmvnorm(t(mom), rep(0, K), M, log = TRUE)  
  alpha <- min(1, exp(logPost_Mom_star-logPost_Mom_t))
  u <- runif(1)
  if(u <= alpha){
    thetaNew <- c(theta)
  }else{
    thetaNew <- thetat
  }
  rest <- list(theta = thetaNew, Prob = alpha)
  return(rest)
}
# Posterior draws HMC
S <- 5000; burnin <- 500; tot <- S + burnin
epsilon <- 0.0025;  L <- 100; M <- diag(2)
# epsilon <- 0.0032;  L <- 100; M <- solve(VarPost)
thetaPostHMC <- matrix(NA, tot, K)
ProbAcceptHMC  <- rep(NA, tot)
thetaHMC <- c(0, 0)
for(s in 1:tot){
  ResHMC <- HMC(theta = thetaHMC, epsilon = epsilon, L = L, M = M)
  thetaHMC <- ResHMC$theta
  thetaPostHMC[s,] <- thetaHMC
  ProbAcceptHMC[s] <- ResHMC$Prob
}
keep <- burnin:S
summary(ProbAcceptHMC[keep])
thetaPostHMCMCMC <- coda::mcmc(thetaPostHMC[keep,])
plot(thetaPostHMCMCMC); coda::autocorr.plot(thetaPostHMCMCMC)
plot(thetaPostHMC[keep,], xlim=range(theta1), ylim=range(theta2), pch=16, col=grey(0.8),
     xlab=expression(theta[1]), ylab=expression(theta[2]))
contour(theta1,theta2,f, drawlabels=FALSE, add=TRUE, col = "blue", lwd = 1.2)
title("Hamiltonian Monte Carlo")

########################## Importance sampling: Beta distribution ########################## 
rm(list = ls()); set.seed(010101); S <- 50000
# Importance sampling from standard normal proposal 
an <- 16.55; bn <- 39.57 # Posterior parameters
theta <- runif(S) # Proposal
ws <- dbeta(theta, an, bn) # Weights
wstars <- ws/sum(ws) # Standardized weights
L <- 20000 # Size of posterior sample
thetaBeta <- sample(theta, L, replace = TRUE, prob = wstars) # Posterior draws
# Figure
h <- hist(thetaBeta, breaks=50, col="blue", xlab="x", main="Beta draws from importance sampling: Uniform (0,1) proposal")
pfit <- seq(min(thetaBeta),max(thetaBeta),length=50)
yfit<-dbeta(pfit, an, bn)
yfit <- yfit*diff(h$mids[1:2])*length(thetaBeta)
lines(pfit, yfit, col="red", lwd=2)
a0 <- 1.44; b0 <- 2.57 # Hyperparameters
s <- 15; n <- 52 # Data
LogMarLik <- lbeta(an, bn)-lbeta(a0, b0) # Marginal likelihood
LogMarLik
# Gelfand-Day method
LikPrior <- function(theta){
  Liki <- theta^s*(1-theta)^(n-s)
  Priori <- dbeta(theta, a0, b0)
  LikPriori <- 1 / Liki * Priori 
  return(LikPriori)
}
LogMarLikGD <- log(1/mean(sapply(thetaBeta, LikPrior)))
LogMarLikGD

########################## SIS and Particle filter: Gordon (1993) example 1 ########################## 
rm(list = ls()); set.seed(010101)
S <- 50000 # Number of particles
sigma_w <- 10^0.5  # State noise
sigma_mu <- 1  # Observation noise
T <- 200 # Sample size
# Simulate true states and observations
theta_true <- numeric(T)
y_obs <- numeric(T)
theta_true[1] <- rnorm(1, mean = 0, sd = sigma_w)  # Initial state
y_obs[1] <- rnorm(1, mean = theta_true[1]^2/20, sd = sigma_mu) 
for (t in 2:T) {
  meant <- 0.5*theta_true[t-1] + 25*(theta_true[t-1]/(1+theta_true[t-1]^2)) + 8*cos(1.2*t)
  theta_true[t] <-  rnorm(1, mean = meant, sd = sigma_w)
  y_obs[t] <- rnorm(1, mean = theta_true[t]^2/20, sd = sigma_mu)  
}
plot(theta_true, type = "l"); plot(y_obs, type = "l")
# Sequential Importance Sampling (SIS)
particles <- matrix(0, nrow = S, ncol = T)  # Store particles
weights <- matrix(0, nrow = S, ncol = T)   # Store weights
weightsSt <- matrix(0, nrow = S, ncol = T)   # Store standardized weights
# Initialization
particles[, 1] <- rnorm(S, mean = 0, sd = sigma_w)  # Sample initial particles
weights[, 1] <- dnorm(y_obs[1], mean = particles[, 1]^2/20, sd = sigma_mu)  # Importance weights
weightsSt[, 1] <- weights[, 1] / sum(weights[, 1])  # Normalize weights
# Sequential updating
for (t in 2:T) {
  # Propagate particles
  particles[, t] <- rnorm(S, mean = 0.5*particles[, t-1] + 25*(particles[, t-1]/(1+particles[, t-1]^2)) + 8*cos(1.2*t), sd = sigma_w)  # Compute weights
  weights[, t] <- weightsSt[, t-1] * dnorm(y_obs[t], mean = particles[, t]^2/20, sd = sigma_mu)
  weightsSt[, t] <- weights[, t] / sum(weights[, t])  # Normalize weights
}
# Estimate the states (weighted mean)
FilterDist <- colSums(particles * weightsSt)
SDFilterDist <- (colSums(particles^2 * weightsSt) - FilterDist^2)^0.5
MargLik <- colMeans(weights)
library(dplyr)
library(ggplot2)
require(latex2exp)
ggplot2::theme_set(theme_bw())
df <- tibble(t = 1:T,
             mean = FilterDist,
             lower = FilterDist - 2*SDFilterDist,
             upper = FilterDist + 2*SDFilterDist,
             theta_true = theta_true)
# Function to plot
plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1,
                fill = "lightblue") +
    geom_line(aes(y = theta_true), colour = "black", alpha = 1,
              linewidth = 0.5) +
    geom_line(aes(y = mean), colour = "blue", linewidth = 0.5) +
    ylab(TeX("$\\theta_{t}$")) + xlab("Time")
  print(p)
}
plot_filtering_estimates(df)

# Particle filtering
particles <- matrix(0, nrow = S, ncol = T)  # Store particles
particlesT <- matrix(0, nrow = S, ncol = T)  # Store resampling particles
weights <- matrix(0, nrow = S, ncol = T)   # Store weights
weightsSt <- matrix(0, nrow = S, ncol = T)   # Store standardized weights
weightsSTT <- matrix(1/S, nrow = S, ncol = T)   # Store standardized weights
logalphas <- matrix(0, nrow = S, ncol = T)   # Store log incremental weights
# Initialization
particles[, 1] <- rnorm(S, mean = 0, sd = sigma_w)  # Sample initial particles
weights[, 1] <- dnorm(y_obs[1], mean = particles[, 1]^2/20, sd = sigma_mu)  # Importance weights
weightsSt[, 1] <- weights[, 1] / sum(weights[, 1])  # Normalize weights
ind <- sample(1:S, size = S, replace = TRUE, prob = weightsSt[, 1]) # Resample 
particles[, 1] <- particles[ind, 1] # Resampled particles
particlesT[, 1] <- particles[, 1] # Resampled particles
# Sequential updating
pb <- winProgressBar(title = "progress bar", min = 0, max = T, width = 300)
for (t in 2:T) {
  # Propagate particles
  particles[, t] <- rnorm(S, mean = 0.5*particles[, t-1] + 25*(particles[, t-1]/(1+particles[, t-1]^2)) + 8*cos(1.2*t), sd = sigma_w)
  # Compute weights
  logalphas[, t] <- dnorm(y_obs[t], mean = particles[, t]^2/20, sd = sigma_mu, log = TRUE) 
  weights[, t] <- exp(logalphas[, t])
  weightsSt[, t] <- weights[, t] / sum(weights[, t])
  if(t < T){
    ind <- sample(1:S, size = S, replace = TRUE, prob = weightsSt[, t])
    particles[, 1:t] <- particles[ind, 1:t]
  }else{
    ind <- sample(1:S, size = S, replace = TRUE, prob = weightsSt[, t])
    particlesT[, 1:t] <- particles[ind, 1:t]
  }
  setWinProgressBar(pb, t, title=paste( round(t/T*100, 0), "% done"))
}
close(pb)
FilterDist <- colSums(particles * weightsSt)
SDFilterDist <- (colSums(particles^2 * weightsSt) - FilterDist^2)^0.5
FilterDistT <- colSums(particlesT * weightsSTT)
SDFilterDistT <- (colSums(particlesT^2 * weightsSTT) - FilterDistT^2)^0.5
library(dplyr)
library(ggplot2)
require(latex2exp)
ggplot2::theme_set(theme_bw())
df <- tibble(t = 1:T,
             mean = FilterDist,
             lower = FilterDist - 2*SDFilterDist,
             upper = FilterDist+ 2*SDFilterDist,
             meanT = FilterDistT,
             lowerT = FilterDistT - 2*SDFilterDistT,
             upperT = FilterDistT + 2*SDFilterDistT,
             x_true = theta_true)

plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1,
                fill = "lightblue") +
    geom_line(aes(y = x_true), colour = "black", alpha = 1,
              linewidth = 0.5) +
    geom_line(aes(y = mean), colour = "blue", linewidth = 0.5) +
    geom_line(aes(y = meanT), colour = "purple", linewidth = 0.5) +
    ylab(TeX("$\\theta_{t}$")) + xlab("Time")
  print(p)
}
plot_filtering_estimates(df)
MargLik <- colMeans(weights)
plot(MargLik, type = "l", xlab = "Time", ylab = "Marginal likelihood", main = "Marginal likelihood")
matplot(1:T, t(weights), type = "l", lty = 1, col = rgb(0, 0, 0, 0.1), main = "Particle Weights Over Time", ylab = "Weight", xlab = "Time")
max(weightsSt[,T])

########################## Metropolis-Hastings: Beta posterior########################## 
rm(list = ls())
set.seed(010101)
an <- 16.55; bn <- 39.57
S <- 100000
p <- runif(S)
accept <- rep(0, S)
for (s in 2:S){
  pc <- runif(1)
  a <- dbeta(pc, an, bn)/dbeta(p[s-1], an, bn)
  U <- runif(1)
  if(U <= a){
    p[s] <- pc
    accept[s] <- 1
  }else{
    p[s] <- p[s-1]
    accept[s] <- 0
  }
}
mean(accept)
mean(p); sd(p)
an/(an + bn); (an*bn/((an+bn)^2*(an+bn+1)))^0.5
h <- hist(p, breaks=50, col="blue", xlab="Proportion Ph.D. students sleeping at least 6 hours", main="Beta draws from a Metropolis-Hastings algorithm")
pfit <- seq(min(p),max(p),length=50)
yfit<-dbeta(pfit, an, bn)
yfit <- yfit*diff(h$mids[1:2])*length(p)
lines(pfit, yfit, col="red", lwd=2)
library(coda); library(latex2exp)
thetaPost <- mcmc(p)
plot(thetaPost, density = FALSE, main = "Trace plot", ylab = "Proportion")
autocorr.plot(thetaPost, main = "Autocorrelation plot")
summary(thetaPost); effectiveSize(thetaPost) 
geweke.diag(thetaPost, frac1 = 0.1, frac2 = 0.5)
raftery.diag(thetaPost, q = 0.025, r = 0.005, s = 0.95)
heidel.diag(thetaPost, eps = 0.1, pvalue = 0.05)
burnin <- 20000; thin <- 10
keep <- seq(burnin, S, thin)
thetaPostNew <- mcmc(p[keep])
plot(thetaPostNew, density = FALSE, main = "Trace plot", ylab = "Proportion")
autocorr.plot(thetaPostNew, main = "Autocorrelation plot")
summary(thetaPostNew); effectiveSize(thetaPostNew) 
geweke.diag(thetaPostNew, frac1 = 0.1, frac2 = 0.5)
raftery.diag(thetaPostNew, q = 0.025, r = 0.005, s = 0.95)
heidel.diag(thetaPostNew, eps = 0.1, pvalue = 0.05)
# Geweke 2004
# Marginal-conditional simulator
a0 <- 1.22; b0 <- 2.57
ThetaPrior <- rbeta(S,a0,b0) 
yPrior <- rbinom(S, 1, prob = ThetaPrior)
parsmcmc1 <- coda::mcmc(cbind(ThetaPrior, ThetaPrior^2))
Summ1 <- summary(parsmcmc1)
# Successive-conditional simulator
N <- 52
SucConSim <- function(a0, b0, par){
  y <- rbinom(N, 1, prob = par)
  an <- a0 + sum(y); bn <- b0 + N - sum(y) 
  pc <- runif(1)
  a <- dbeta(pc, an, bn)/dbeta(par, an, bn)
  U <- runif(1)
  if(U <= a){
    par <- pc
  }else{
    par <- par
  }
  return(par)
}
# a0 <- 1.22; b0 <- 2.57
a0 <- 2.57; b0 <- 1.22
par2 <- rbeta(1, a0, b0)
pars2 <- rep(NA, S)
pars2[1] <- par2
for(s in 2:S){
  Res <- SucConSim(a0 = a0, b0 = b0, par = pars2[s-1])
  pars2[s] <- Res
}
parsmcmc2 <- coda::mcmc(cbind(pars2, pars2^2))
Summ2 <- summary(parsmcmc2)
TestGeweke <- function(j){
  Test <- (Summ1[["statistics"]][j,1] - Summ2[["statistics"]][j,1])/(Summ1[["statistics"]][j,4]+Summ2[["statistics"]][j,4])^0.5
  Reject <- abs(Test) > qnorm(0.975)
  return(list(Test = Test, Reject = Reject))
}
TestGeweke(1); TestGeweke(2)
