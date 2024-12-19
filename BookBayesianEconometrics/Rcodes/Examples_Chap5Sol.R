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
