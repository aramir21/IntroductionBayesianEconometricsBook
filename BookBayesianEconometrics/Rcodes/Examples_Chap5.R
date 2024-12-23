########################## Mining change point ########################## 
rm(list = ls())
set.seed(010101)
dataset<-read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/MiningDataCarlin.csv",header=T)
attach(dataset)
str(dataset)
a10 <- 0.5; a20 <- 0.5
b10 <- 1; b20 <- 1
y <- Count
sumy <- sum(Count); T <- length(Count)
theta1 <- NULL; theta2 <- NULL
kk <- NULL; H <- 60
MCMC <- 20000; burnin <- 1000; S <- MCMC + burnin; keep <- (burnin+1):S
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
for(s in 1:S){
  a1 <- a10 + sum(y[1:H])
  b1 <- b10+H
  theta11 <- rgamma(1,a1,b1)
  theta1 <- c(theta1,theta11)
  a2 <- a20 + sum(y[(1+H):T])
  b2 <- b20 + T-H
  theta22 <- rgamma(1,a2,b2)
  theta2 <- c(theta2,theta22)
  pp<-NULL
  for(l in 1:T){
    p <- exp(l*(theta22-theta11))*(theta11/theta22)^(sum(y[1:l]))
    pp <- c(pp,p)
  }
  prob <- pp/sum(pp)
  H <- sample(1:T,1,prob=prob)
  kk <- c(kk,H)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
library(coda)
theta1Post <- mcmc(theta1[keep])
summary(theta1Post)
theta2Post <- mcmc(theta2[keep])
summary(theta2Post)
HPost <- mcmc(kk); summary(HPost)
hist(HPost, main = "Histogram: Posterior mean change point", xlab = "Posterior mean", col = "blue", breaks = 25)

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

########################## Hamiltonian Monte Carlo: Bi-variate normal distribution ########################## 
rm(list = ls()); set.seed(010101)
# Gibbs sampler
Gibbs <- function(theta, rho){
  thetal <- rnorm(1, mean = rho*theta, sd = (1- rho^2)^0.5)
  return(thetal)
}
# Metropolis-Hastings
MH <- function(theta, rho, sig2){
  SIGMA <- matrix(c(1, rho, rho, 1), 2, 2)
  SIGMAc <- matrix(c(1, sig2, sig2, 1), 2, 2)
  thetac <- MASS::mvrnorm(1, mu = theta, Sigma = SIGMAc)
  a <- mvtnorm::dmvnorm(thetac, c(0, 0), SIGMA)/mvtnorm::dmvnorm(theta, c(0, 0), SIGMA)
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
# Hamiltonian Monte Carlo
HMC <- function(theta, rho, epsilon, M){
  SIGMA <- matrix(c(1, rho, rho, 1), 2, 2) 
  L <- ceiling(1/epsilon)
  Minv <- solve(M); thetat <- theta
  K <- length(thetat)
  mom <- t(mvtnorm::rmvnorm(1, rep(0, K), M))
  logPost_Mom_t <- mvtnorm::dmvnorm(t(theta), rep(0, K), SIGMA, log = TRUE) +  mvtnorm::dmvnorm(t(mom), rep(0, K), M, log = TRUE)  
  for(l in 1:L){
    if(l == 1 | l == L){
      mom <- mom + 0.5*epsilon*(-solve(SIGMA)%*%theta)
      theta <- theta + epsilon*Minv%*%mom
    }else{
      mom <- mom + epsilon*(-solve(SIGMA)%*%theta)
      theta <- theta + epsilon*Minv%*%mom
    }
  }
  logPost_Mom_star <- mvtnorm::dmvnorm(t(theta), rep(0, K), SIGMA, log = TRUE) +  mvtnorm::dmvnorm(t(mom), rep(0, K), M, log = TRUE)  
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
# Hyperparameters
rho <- 0.98; sig2 <- 0.18^2
# Posterior draws Gibbs and M-H
S <- 8000; thin <- 20; K <- 2
thetaPostGibbs <- matrix(NA, S, K)
thetaPostMH <- matrix(NA, S, K)
AcceptMH <- rep(NA, S)
thetaGibbs <- c(-2, 3); thetaMH <- c(-2, 3)
for(s in 1:S){
  theta1 <- Gibbs(thetaGibbs[2], rho)
  theta2 <- Gibbs(theta1, rho)
  thetaGibbs <- c(theta1, theta2)
  ResMH <- MH(thetaMH, rho, sig2)
  thetaMH <- ResMH$theta
  thetaPostGibbs[s,] <- thetaGibbs
  thetaPostMH[s,] <- thetaMH
  AcceptMH[s] <- ResMH$accept
}
keep <- seq(0, S, thin)
mean(AcceptMH[keep[-1]])
thetaPostGibbsMCMC <- coda::mcmc(thetaPostGibbs[keep[-1],])
summary(thetaPostGibbsMCMC)
coda::autocorr.plot(thetaPostGibbsMCMC)
thetaPostMHMCMC <- coda::mcmc(thetaPostMH[keep[-1],])
plot(thetaPostMHMCMC)
coda::autocorr.plot(thetaPostMHMCMC)
# Posterior draws HMC
S <- 400;epsilon <- 0.05;  L <- ceiling(1/epsilon); M <- diag(2)
thetaPostHMC <- matrix(NA, S, K)
ProbAcceptHMC  <- rep(NA, S)
thetaHMC <- c(-2, 3)
for(s in 1:S){
  ResHMC <- HMC(theta = thetaHMC, rho, epsilon, M)
  thetaHMC <- ResHMC$theta
  thetaPostHMC[s,] <- thetaHMC
  ProbAcceptHMC[s] <- ResHMC$Prob
}
thetaPostHMCMCMC <- coda::mcmc(thetaPostHMC)
plot(thetaPostHMCMCMC); coda::autocorr.plot(thetaPostHMCMCMC)
summary(ProbAcceptHMC)

sd(thetaPostHMC[,1]); sd(thetaPostMH[keep[-1],1]); sd(thetaPostGibbs[keep[-1],1])
#Figure
df <- as.data.frame(cbind(1:S, thetaPostHMC[,1], thetaPostMH[keep[-1],1], thetaPostGibbs[keep[-1],1]))
colnames(df) <- c("Iter", "HMC", "MH", "Gibbs")
require(latex2exp)
g1 <- ggplot(df, aes(x= Iter)) + geom_point(aes(y=HMC), colour="black") + labs(x = "Iteration", y = TeX("$\\theta_{1}$"), title = "HMC algorithm")
g2 <- ggplot(df, aes(x= Iter)) + geom_point(aes(y=MH), colour="black") + labs(x = "Iteration", y = TeX("$\\theta_{1}$"), title = "M-H algorithm")
g3 <- ggplot(df, aes(x= Iter)) + geom_point(aes(y=Gibbs), colour="black") + labs(x = "Iteration", y = TeX("$\\theta_{1}$"), title = "Gibbs sampling")
library(ggpubr)
ggarrange(g3, g2, g1, labels = c("A", "B", "C"), ncol = 3, nrow = 1)

########################## Importance sampling: Cauchy distribution ########################## 
rm(list = ls()); set.seed(010101)
S <- 20000 # Size proposal
# Importance sampling from standard normal proposal 
thetaNs <- rnorm(S)
wNs <- dcauchy(thetaNs)/dnorm(thetaNs)
wNstars <- wNs/sum(wNs)
L <- 10000 # Size posterior
thetaCauchyN <- sample(thetaNs, L, replace = TRUE, prob = wNstars)
h <- hist(thetaCauchyN, breaks=50, col="blue", xlab="x", main="Cauchy draws from importance sampling: Normal standard proposal")
pfit <- seq(min(thetaCauchyN),max(thetaCauchyN),length=50)
yfit<-dcauchy(pfit)
yfit <- yfit*diff(h$mids[1:2])*length(thetaCauchyN)
lines(pfit, yfit, col="red", lwd=2)
# Importance sampling from Student's t proposal 
df <- 3
thetaTs <- rt(S, df = df)
wTs <- dcauchy(thetaTs)/dt(thetaTs, df = df)
wTstars <- wTs/sum(wTs)
thetaCauchyT <- sample(thetaTs, L, replace = TRUE, prob = wTstars)
h <- hist(thetaCauchyT, breaks=50, col="blue", xlab="x", main="Cauchy draws from importance sampling: Student's t proposal")
pfit <- seq(min(thetaCauchyT),max(thetaCauchyT),length=50)
yfit<-dcauchy(pfit)
yfit <- yfit*diff(h$mids[1:2])*length(thetaCauchyT)
lines(pfit, yfit, col="red", lwd=2)
plot(wNstars, main = "Importance sampling: Cauchy distribution", ylab = "Weights", xlab = "Iterations")
points(wTstars, col = "blue")
legend("topright", legend = c("Normal", "Student's t"), col = c("black", "blue"), pch = c(1, 1))

########################## Sequential importance sampling: Dynamic linear model ########################## 
rm(list = ls()); set.seed(010101)
S <- 50000 # Number of particles
sigma_w <- 1  # State noise
sigma_mu <- 1  # Observation noise
phi <- 0.5 # Coefficient in observation equation
T <- 50 # Sample size
# Simulate true states and observations
theta_true <- numeric(T)
y_obs <- numeric(T)
theta_true[1] <- rnorm(1, mean = 0, sd = sigma_w)  # Initial state
for (t in 2:T) {
  theta_true[t] <- rnorm(1, mean = theta_true[t-1], sd = sigma_w)
}
y_obs <- rnorm(T, mean = phi*theta_true, sd = sigma_mu)
# Sequential Importance Sampling (SIS)
particles <- matrix(0, nrow = S, ncol = T)  # Store particles
weights <- matrix(0, nrow = S, ncol = T)   # Store weights
weightsSt <- matrix(0, nrow = S, ncol = T)   # Store standardized weights
# Initialization
particles[, 1] <- rnorm(S, mean = 0, sd = sigma_w)  # Sample initial particles
weights[, 1] <- dnorm(y_obs[1], mean = phi*particles[, 1], sd = sigma_mu)  # Importance weights
weightsSt[, 1] <- weights[, 1] / sum(weights[, 1])  # Normalize weights
# Sequential updating
for (t in 2:T) {
  # Propagate particles
  particles[, t] <- rnorm(S, mean = particles[, t-1], sd = sigma_w)
  # Compute weights
  weights[, t] <- weightsSt[, t-1] * dnorm(y_obs[t], mean = phi*particles[, t], sd = sigma_mu) # Recursive weight update
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

matplot(1:T, t(weights), type = "l", lty = 1, col = rgb(0, 0, 0, 0.1), main = "Particle Weights Over Time", ylab = "Weight", xlab = "Time")

########################## Particle filter: Dynamic linear model ########################## 
rm(list = ls()); set.seed(010101)
S <- 50000 # Number of particles
sigma_w <- 1  # State noise
sigma_mu <- 1  # Observation noise
phi <- 0.5 # Coefficient in observation equation
T <- 200 # Sample size
# Simulate true states and observations
theta_true <- numeric(T)
y_obs <- numeric(T)
theta_true[1] <- rnorm(1, mean = 0, sd = sigma_w)  # Initial state
for (t in 2:T) {
  theta_true[t] <- rnorm(1, mean = theta_true[t-1], sd = sigma_w)
}
y_obs <- rnorm(T, mean = phi*theta_true, sd = sigma_mu)
# Particle filtering
particles <- matrix(0, nrow = S, ncol = T)  # Store particles
particlesT <- matrix(0, nrow = S, ncol = T)  # Store resampling particles
weights <- matrix(0, nrow = S, ncol = T)   # Store weights
weightsSt <- matrix(0, nrow = S, ncol = T)   # Store standardized weights
weightsSTT <- matrix(1/S, nrow = S, ncol = T)   # Store standardized weights
logalphas <- matrix(0, nrow = S, ncol = T)   # Store log incremental weights
# Initialization
particles[, 1] <- rnorm(S, mean = 0, sd = sigma_w)  # Sample initial particles
weights[, 1] <- dnorm(y_obs[1], mean = phi*particles[, 1], sd = sigma_mu)  # Importance weights
weightsSt[, 1] <- weights[, 1] / sum(weights[, 1])  # Normalize weights
ind <- sample(1:S, size = S, replace = TRUE, prob = weightsSt[, 1]) # Resample 
particles[, 1] <- particles[ind, 1] # Resampled particles
particlesT[, 1] <- particles[, 1] # Resampled particles
# Sequential updating
pb <- winProgressBar(title = "progress bar", min = 0, max = T, width = 300)
for (t in 2:T) {
  # Propagate particles
  particles[, t] <- rnorm(S, mean = particles[, t-1], sd = sigma_w)
  # Compute weights
  logalphas[, t] <- dnorm(y_obs[t], mean = phi*particles[, t], sd = sigma_mu, log = TRUE) 
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
plot(MargLik, type = "l")
matplot(1:T, t(weights), type = "l", lty = 1, col = rgb(0, 0, 0, 0.1), main = "Particle Weights Over Time", ylab = "Weight", xlab = "Time")

########################## Mining change point ########################## 
rm(list = ls())
set.seed(010101)
dataset<-read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/MiningDataCarlin.csv",header=T)
attach(dataset)
str(dataset)
a10 <- 0.5; a20 <- 0.5
b10 <- 1; b20 <- 1
y <- Count
sumy <- sum(Count); T <- length(Count)
theta1 <- NULL; theta2 <- NULL
kk <- NULL; H <- 60
MCMC <- 20000; burnin <- 1000; S <- MCMC + burnin; keep <- (burnin+1):S
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
for(s in 1:S){
  a1 <- a10 + sum(y[1:H])
  b1 <- b10+H
  theta11 <- rgamma(1,a1,b1)
  theta1 <- c(theta1,theta11)
  a2 <- a20 + sum(y[(1+H):T])
  b2 <- b20 + T-H
  theta22 <- rgamma(1,a2,b2)
  theta2 <- c(theta2,theta22)
  pp<-NULL
  for(l in 1:T){
    p <- exp(l*(theta22-theta11))*(theta11/theta22)^(sum(y[1:l]))
    pp <- c(pp,p)
  }
  prob <- pp/sum(pp)
  H <- sample(1:T,1,prob=prob)
  kk <- c(kk,H)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
library(coda); library(latex2exp)
theta1Post <- mcmc(theta1[keep]); summary(theta1Post)
HPost <- mcmc(kk); summary(HPost)
theta2Post <- mcmc(theta2[keep])
plot(theta2Post, density = FALSE, main = "Trace plot", ylab = TeX("$\\lambda_{2}$"))
autocorr.plot(theta2Post, main = "Autocorrelation plot")
summary(theta2Post); effectiveSize(theta2Post) 
geweke.diag(theta2Post, frac1 = 0.1, frac2 = 0.5)
raftery.diag(theta2Post, q = 0.025, r = 0.005, s = 0.95)
heidel.diag(theta2Post, eps = 0.1, pvalue = 0.05)


