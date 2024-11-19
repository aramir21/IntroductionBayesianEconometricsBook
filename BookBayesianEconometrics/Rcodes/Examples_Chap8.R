########################## Simulation exercise: Dynamic linear model ########################## 
rm(list = ls()); set.seed(010101)
T <- 200; sig2 <- 0.5^2
x <- rnorm(T, mean = 1, sd = 1) 
X <- cbind(1, x); B0 <- c(1, 0.5)
K <- length(B0)
e <- rnorm(T, mean = 0, sd = sig2^0.5)
Omega <- diag(c(0.2, 0.1))
w <- MASS::mvrnorm(T, c(0, 0), Omega)
Bt <- matrix(NA, T, K); Bt[1,] <- B0
yt <- rep(NA, T) 
yt[1] <- X[1,]%*%B0 + e[1]
yt[1] <- x[1]*Bt[1] + e[1]
for(t in 1:T){
  if(t == 1){
    Bt[t,] <- w[t,]
  }else{
    Bt[t,] <- Bt[t-1,] + w[t,]
  }
  yt[t] <- X[t,]%*%Bt[t,] + e[t]
}
RegLS <- lm(yt ~ x)
SumRegLS <- summary(RegLS)
SumRegLS; SumRegLS$sigma^2  
Bp <- matrix(RegLS$coefficients, T, K, byrow = TRUE)
S <- 20
for(t in S:T){
  RegLSt <- lm(yt[1:t] ~ x[1:t])
  Bp[t,] <- RegLSt$coefficients 
}
# plot(Bp[S:T,2], type = "l")
VarBp <- var(Bp)
# State spece model
ModelReg <- function(par){
  Mod <- dlm::dlmModReg(x, dV = exp(par[1]), dW = exp(par[2:3]), m0 = RegLS$coefficients,
                        C0 = VarBp)
  return(Mod)
}
outMLEReg <- dlm::dlmMLE(yt, parm = rep(0, K+1), ModelReg)
exp(outMLEReg$par)
RegFilter <- dlm::dlmFilter(yt, ModelReg(outMLEReg$par))
RegSmoth <- dlm::dlmSmooth(yt, ModelReg(outMLEReg$par))
SmoothB2 <- RegSmoth$s[-1,2]
VarSmooth <- dlm::dlmSvd2var(u = RegSmoth[["U.S"]], RegSmoth[["D.S"]])
SDVarSmoothB2 <- sapply(2:(T+1), function(t){VarSmooth[[t]][K,K]^0.5}) 
LimInfB2 <- SmoothB2 - qnorm(0.975)*SDVarSmoothB2
LimSupB2 <- SmoothB2 + qnorm(0.975)*SDVarSmoothB2
# Gibbs
MCMC <- 2000; burnin <- 1000
gibbsOut <- dlm::dlmGibbsDIG(yt, mod = dlm::dlmModReg(x), a.y = SumRegLS$sigma^2, b.y = 10*SumRegLS$sigma^2, a.theta = max(diag(VarBp)), b.theta = 10*max(diag(VarBp)), n.sample = MCMC, thin = 5, save.states = TRUE)
B2t <- matrix(0, MCMC - burnin, T + 1)
for(t in 1:(T+1)){
  B2t[,t] <- gibbsOut[["theta"]][t,2,-c(1:burnin)] 
}
Lims <- apply(B2t, 2, function(x){quantile(x, c(0.025, 0.975))})
summary(coda::mcmc(gibbsOut[["dV"]]))
summary(coda::mcmc(gibbsOut[["dW"]]))
# Figure
require(latex2exp) # LaTeX equations in figures
xx <- c(1:(T+1), (T+1):1)
yy <- c(Lims[1,], rev(Lims[2,]))
plot   (xx, yy, type = "n", xlab = "Time", ylab = TeX("$\\beta_{t2}$"))
polygon(xx, yy, col = "lightblue", border = "lightblue")
xxML <- c(1:T, T:1)
yyML <- c(LimInfB2, rev(LimSupB2))
polygon(xxML, yyML, col = "blue", border = "blue")
lines(colMeans(B2t), col = "red", lw = 2)
lines(Bt[,2], col = "black", lw = 2)
lines(SmoothB2, col = "green", lw = 2)
title("State vector: Slope parameter")

########################## Application of dynamic linear models: Effects of inflation on interest rate ########################## 
rm(list = ls())
set.seed(010101)
DataIntRate <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/16INTDEF.csv", sep = ",", header = TRUE, quote = "")
attach(DataIntRate)
Xt <- cbind(diff(inf), diff(def))
K <- dim(Xt)[2] + 1
yt <- diff(i3)
T <- length(yt)
# plot(yt, type = "l")
# plot(Xt[,2], type = "l")
RegLS <- lm(yt ~ Xt)
SumRegLS <- summary(RegLS)
SumRegLS; SumRegLS$sigma^2  
# Recursive OLS
Bp <- matrix(RegLS$coefficients, T, K, byrow = TRUE)
S <- 20
for(t in S:T){
  RegLSt <- lm(yt[1:t] ~ Xt[1:t,])
  Bp[t,] <- RegLSt$coefficients 
}
# plot(Bp[S:T,2], type = "l")
VarBp <- var(Bp)
# State spece model
ModelReg <- function(par){
  Mod <- dlm::dlmModReg(Xt, dV = exp(par[1]), dW = exp(par[2:(K+1)]), m0 = RegLS$coefficients,
                        C0 = diag(VarBp))
  return(Mod)
}
MCMC <- 12000; burnin <- 2000; thin <- 10
gibbsOut <- dlm::dlmGibbsDIG(yt, mod = dlm::dlmModReg(Xt), a.y = SumRegLS$sigma^2, b.y = 10*SumRegLS$sigma^2, a.theta = max(diag(VarBp)), b.theta = 10*max(diag(VarBp)), n.sample = MCMC, thin = 5, save.states = TRUE)
B2t <- matrix(0, MCMC - burnin, T + 1)
for(t in 1:(T+1)){
  B2t[,t] <- gibbsOut[["theta"]][t,2,-c(1:burnin)] 
}
Lims <- apply(B2t, 2, function(x){quantile(x, c(0.025, 0.975))})
dV <- coda::mcmc(gibbsOut[["dV"]][-c(1:burnin)])
dW <- coda::mcmc(gibbsOut[["dW"]][-c(1:burnin),])
summary(dV)
summary(dW)
plot(dV)
plot(dW)
library(fanplot)
require(latex2exp)
df <- as.data.frame(B2t)
plot(NULL, main="Percentiles", xlim = c(1, T+1), ylim = c(-1, 2), xlab = "Time", ylab = TeX("$\\beta_{t2}$"))
fan(data = df)
lines(colMeans(B2t), col = "black", lw = 2)
abline(h=0, col = "blue")

########################## Application of linear regression with AR errors: Effects of inflation on interest rate ########################## 
rm(list = ls())
set.seed(010101)
DataIntRate <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/16INTDEF.csv", sep = ",", header = TRUE, quote = "")
attach(DataIntRate)
yt <- diff(i3); ytlag <- dplyr::lag(yt, n = 1)
T <- length(yt)
Xt <- cbind(diff(inf), diff(def)); Xtlag <- dplyr::lag(Xt, n = 1)
K <- dim(Xt)[2] + 1

######### Simulation ########
# T <- 500
# B <- c(0.6, -0.7, 1.2); sig2 <- 0.8; Phi <- 0.4
# Xt <- cbind(rnorm(T), rnorm(T)); Xtlag <- dplyr::lag(Xt, n = 1)
# e <- rnorm(T, 0, sig2^0.5)
# mu <- rnorm(T, 0, sd = (sig2/(1-Phi))^0.5)
# for(t in 2:T){
#   mu[t] <- Phi*mu[t-1] + e[t]
# }
# Xt <- cbind(1, Xt)
# yt <- Xt%*%B + mu; ytlag <- dplyr::lag(yt, n = 1)
Reg <- lm(yt ~ ytlag + I(Xt[,-1] - Xtlag))
SumReg <- summary(Reg); SumReg
# Xt <- Xt[,-1]
PostSig2 <- function(Beta, Phi){
  Xstar<- matrix(NA, T-1, K - 1)
  ystar <- matrix(NA, T-1, 1)
  for(t in 2:T){
    Xstar[t-1,] <- Xt[t,] - Phi*Xt[t-1,]
    ystar[t-1,] <- yt[t] - Phi*yt[t-1]
  }
  Xstar <- cbind(1, Xstar)
  an <- T - 1 + a0
  dn <- d0 + t(ystar - Xstar%*%Beta)%*%(ystar - Xstar%*%Beta) + t(Beta - b0)%*%B0i%*%(Beta - b0)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBeta <- function(sig2, Phi){
  Xstar<- matrix(NA, T-1, K - 1)
  ystar <- matrix(NA, T-1, 1)
  for(t in 2:T){
    Xstar[t-1,] <- Xt[t,] - Phi*Xt[t-1,]
    ystar[t-1,] <- yt[t] - Phi*yt[t-1]
  }
  Xstar <- cbind(1, Xstar)
  XtXstar <- t(Xstar)%*%Xstar
  Xtystar <- t(Xstar)%*%ystar
  Bn <- solve(B0i + XtXstar)
  bn <- Bn%*%(B0i%*%b0 + Xtystar)
  # solve(XtXstar)%*%Xtystar
  Beta <- MASS::mvrnorm(1, bn, sig2*Bn)
  return(Beta)
}
PostPhi <- function(sig2, Beta){
  u <- yt - cbind(1,Xt)%*%Beta
  U <- u[-T]
  ustar <- u[-1]
  UtU <- t(U)%*%U
  Utu <- t(U)%*%ustar
  Phin <- solve(Phi0i + sig2^(-1)*UtU)
  phin <- Phin%*%(Phi0i%*%phi0 + sig2^(-1)*Utu)
  Phi <- truncnorm::rtruncnorm(1, a = -1, b = 1, mean = phin, sd = Phin^0.5)
  return(Phi)
}
# Hyperparameters
d0 <- 0.01; a0 <- 0.01
b0 <- rep(0, K); c0 <- 1; 
B0 <- c0*diag(K); B0i <- solve(B0)
phi0 <- 0; Phi0 <- 1; Phi0i <- 1/Phi0
# MCMC parameters
mcmc <- 15000
burnin <- 5000
tot <- mcmc + burnin
thin <- 1
PostBetas <- matrix(0, mcmc+burnin, K)
PostSigma2s <- rep(0, mcmc+burnin)
PostPhis <- rep(0, mcmc+burnin)
Beta <- rep(0, K); Phi <- 0
sig2 <- SumReg$sigma^2; Phi <- SumReg$coefficients[2,1]
Beta <- SumReg$coefficients[c(1,3,4),1]
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  sig2 <- PostSig2(Beta = Beta, Phi = Phi)
  PostSigma2s[s] <- sig2
  Beta <- PostBeta(sig2 = sig2, Phi = Phi)
  PostBetas[s,] <- Beta
  Phi <- PostPhi(sig2 = sig2, Beta = Beta)
  PostPhis[s] <- Phi
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0), "% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- coda::mcmc(PostBetas[keep,])
summary(PosteriorBetas)
PosteriorSigma2 <- coda::mcmc(PostSigma2s[keep])
summary(PosteriorSigma2)
PosteriorPhi <- coda::mcmc(PostPhis[keep])
summary(PosteriorPhi)
dfBinf <- as.data.frame(PosteriorBetas[,2])
# Basic density
p <- ggplot(dfBinf, aes(x=var1)) + 
  geom_density(color="darkblue", fill="lightblue") +
  geom_vline(aes(xintercept=mean(var1)), color="blue", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=quantile(var1, 0.025)), color="red", linetype="dashed", linewidth=1) +
  geom_vline(aes(xintercept=quantile(var1, 0.975)), color="red", linetype="dashed", linewidth=1) +
  labs(title="Density effect of inflation on interest rate", x="Effect of inflation", y = "Density")
p

########################## Simulation exercise: AR(2) model ########################## 
rm(list = ls())
set.seed(010101)
T <- 200
mu <- 0.5 # 0.007; 
phi1 <- 0.5; phi2 <- 0.3; sig <- 0.5 # 0.035
Ey <- mu/(1-phi1-phi2); Sigy <- sig*((1-phi2)/(1-phi2-phi1^2-phi2*phi1^2-phi2^2+phi2^3))^0.5 
y <- rnorm(T, mean = Ey, sd = Sigy)
e <- rnorm(T, mean = 0, sd = sig)
for(t in 3:T){
  y[t] <- mu + phi1*y[t-1] + phi2*y[t-2] + e[t]
}
mean(y); sd(y)
y <- ts(y, start=c(1820, 1), frequency=1)
plot(y)
iter <- 10000; burnin <- 5000; thin <- 1; tot <- iter + burnin
library(bayesforecast)
sf1 <- bayesforecast::stan_sarima(y, order = c(2, 0, 0), prior_mu0 = normal(0, 1),
                                  prior_ar = normal(0, 1), prior_sigma0 = inverse.gamma(0.01/2, 0.01/2),
                                  seasonal = c(0, 0, 0), iter = tot, warmup = burnin, chains = 1)
keep <- seq(burnin+1, tot, thin)
Postmu <- sf1[["stanfit"]]@sim[["samples"]][[1]][["mu0"]][keep]
Postsig <- sf1[["stanfit"]]@sim[["samples"]][[1]][["sigma0"]][keep]
Postphi1 <- sf1[["stanfit"]]@sim[["samples"]][[1]][["ar0[1]"]][keep]
Postphi2 <- sf1[["stanfit"]]@sim[["samples"]][[1]][["ar0[2]"]][keep]
Postdraws <- coda::mcmc(cbind(Postmu, Postsig, Postphi1, Postphi2))
summary(Postdraws)
# par(mar=c(1,1,1,1))
plot(Postdraws)


########################## Hamiltonian Monte Carlo: AR(2) model ########################## 
rm(list = ls())
set.seed(010101)
T <- 1000; K <- 4
mu <- 0.5 
phi1 <- 0.5; phi2 <- 0.3; sig <- 0.5 
Ey <- mu/(1-phi1-phi2); Sigy <- sig*((1-phi2)/(1-phi2-phi1^2-phi2*phi1^2-phi2^2+phi2^3))^0.5 
y <- rnorm(T, mean = Ey, sd = Sigy)
e <- rnorm(T, mean = 0, sd = sig)
for(t in 3:T){
  y[t] <- mu + phi1*y[t-1] + phi2*y[t-2] + e[t]
}
# Hyperparameters
d0 <- 0.01; a0 <- 0.01
mu0 <- 0; MU0 <- 1
phi0 <- c(0, 0); Phi0 <- diag(2)

LogPost <- function(theta, y){
  mu <- theta[1]
  phi1 <- theta[2]
  phi2 <- theta[3]
  tau <- theta[4]
  sig2 <- exp(tau)
  logLik <- NULL
  for(t in 3:T){
    logLikt <- dnorm(y[t], mean = mu + phi1*y[t-1] + phi2*y[t-2], sd = sig2^0.5, log = TRUE)
    logLik <- c(logLik, logLikt)
  }
  logLik <- sum(logLik)
  # logPrior <- dnorm(mu, mean = mu0, sd = MU0^0.5, log = TRUE) + dnorm(phi1, mean = phi0[1], sd = Phi0[1,1]^0.5, log = TRUE) +
  #   + dnorm(phi2, mean = phi0[2], sd = Phi0[2,2]^0.5, log = TRUE) + a0/2*log(d0/2)-lgamma(a0/2)-a0/2*tau-d0/2*exp(tau)
  # logPosterior <- logLik + logPrior
  logPrior <- dnorm(mu, mean = mu0, sd = MU0^0.5, log = TRUE) + dnorm(phi1, mean = phi0[1], sd = Phi0[1,1]^0.5, log = TRUE) +
    + dnorm(phi2, mean = phi0[2], sd = Phi0[2,2]^0.5, log = TRUE) + invgamma::dinvgamma(sig2, shape = a0/2, rate = d0/2, log = TRUE)
  logPosterior <- logLik + logPrior + tau
  # return(logPosterior)
  return(-logPosterior) # Multiply by -1 to minimize using optim
  # return(-logLik)
}
theta0 <- c(mean(y), 0, 0, var(y))
LogPost(theta = theta0, y)
Opt <- optim(theta0, LogPost, y = y, hessian = TRUE)
theta0 <- Opt$par
VarPost <- solve(Opt$hessian)
# Gradient log posterior
GradientThetaNum <- function(theta, y){
  e <- 0.0001
  GradNum <- NULL
  for(k in 1:length(theta)){
    thetalow <- theta
    thetahigh <- theta
    thetalow[k] <- theta[k] - e
    thetahigh[k] <- theta[k] + e
    LogPostlow <- LogPost(theta = thetalow, y)
    LogPosthigh <- LogPost(theta = thetahigh, y)
    GradNumk <- (LogPosthigh - LogPostlow)/(2*e)
    GradNum <- c(GradNum, GradNumk)
  }
  return(GradNum)
}

GradientTheta <- function(theta, y){
  mu <- theta[1]
  phi1 <- theta[2]
  phi2 <- theta[3]
  tau <- theta[4]
  sig2 <- exp(tau)
  SumLik <- matrix(0, 3, 1)
  SumLik2 <- NULL
  for(t in 3:T){
    xt <- matrix(c(1, y[t-1], y[t-2]), 3, 1)
    SumLikt <- (y[t] - (mu + phi1*y[t-1] + phi2*y[t-2]))*xt
    SumLik2t <- (y[t] - (mu + phi1*y[t-1] + phi2*y[t-2]))^2
    SumLik <- rowSums(cbind(SumLik, SumLikt))
    SumLik2 <- sum(SumLik2, SumLik2t)
  }
  Grad_mu <- SumLik[1]/sig2 - (1/MU0)*(mu - mu0)
  Grad_phi1 <- SumLik[2]/exp(tau) - 1/Phi0[1,1]*(phi1 - phi0[1])
  Grad_phi2 <- SumLik[3]/exp(tau) - 1/Phi0[2,2]*(phi2 - phi0[2])
  Grad_tau <- -(T-2)/2 + SumLik2/(2*exp(tau)) - (a0/2 + 1) + d0/(2*exp(tau)) + 1 
  Grad <- c(Grad_mu, Grad_phi1, Grad_phi2, Grad_tau)
  return(Grad)
}
GradientTheta(theta = Opt$par, y)
StatRest <- function(phi1, phi2){
  if(abs(phi2) < 1 & phi1 + phi2 < 1 & phi2 - phi1 < 1){
    check <- 1
  }else{
    check <- 0
  }
  return(check)
}
StatRest(phi1 = -1.5, phi2 = -0.3)
# Hamiltonian Monte Carlo function
HMC <- function(theta, y, epsilon, M){
  L <- ceiling(1/epsilon)
  Minv <- solve(M); thetat <- theta
  K <- length(thetat)
  mom <- t(mvtnorm::rmvnorm(1, rep(0, K), M))
  logPost_Mom_t <- -LogPost(thetat, y) +  mvtnorm::dmvnorm(t(mom), rep(0, K), M, log = TRUE)  
  for(l in 1:L){
    if(l == L){
      mom <- mom + 0.5*epsilon*GradientTheta(theta, y)
      theta <- theta + epsilon*Minv%*%mom
    }else{
      mom <- mom + epsilon*GradientTheta(theta, y)
      theta <- theta + epsilon*Minv%*%mom
    }
  }
  logPost_Mom_star <- -LogPost(theta, y) +  mvtnorm::dmvnorm(t(mom), rep(0, K), M, log = TRUE)  
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
epsilon <- 0.1; M <- solve(VarPost)
theta <- theta0
HMC(theta, y, epsilon, M)

# Posterior draws
S <- 1000; burnin <- 1000; thin <- 2; tot <- S + burnin
thetaPost <- matrix(NA, tot, K)
ProbAcept <- rep(NA, tot)
# theta0 <- theta0
theta0 <- c(mean(y), 0, 0, exp(var(y)))  
M <- solve(VarPost); epsilon0 <- 0.1
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  epsilon <- runif(1, 0, 2*epsilon0)
  L <- ceiling(1/epsilon)
  HMCs <- HMC(theta = theta0, y, epsilon, M) 
  theta0 <- HMCs$theta 
  thetaPost[s,] <- HMCs$theta
  ProbAcept[s] <- HMCs$Prob
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0), "% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
thetaF <- coda::mcmc(thetaPost[keep,])
summary(thetaF)
plot(thetaF)
summary(exp(thetaF[,K]))
plot(exp(thetaF[,K]))
ProbAceptF <- coda::mcmc(ProbAcept[keep])
summary(ProbAceptF)

########################## Stochastic volatility models: MCMC approach ########################## 
rm(list = ls())
set.seed(010101)
T <- 1250; K <- 2
X <- matrix(rnorm(T*K), T, K)
B <- c(0.5, 0.3); mu <- -10; phi <- 0.95; sigma <- 0.3
h <- numeric(T)
y <- numeric(T)
h[1] <- rnorm(1, mu, sigma / sqrt(1 - phi^2))  # Initial state
y[1] <- X[1,]%*%B + rnorm(1, 0, exp(h[1] / 2))           # Initial observation
for (t in 2:T) {
  h[t] <- mu + phi*(h[t-1]-mu) + rnorm(1, 0, sigma)
  y[t] <- X[t,]%*%B + rnorm(1, 0, sd = exp(0.5*h[t]))
}
plot(y, type = "l") 
plot(h, type = "l")
df <- as.data.frame(cbind(y, X))
colnames(df) <- c("y", "x1", "x2")
MCMC <- 10000; burnin <- 10000; thin <- 5
res <- stochvol::svsample(y, designmatrix = X, draws = MCMC, burnin = burnin, thin = thin, priormu = c(0, 100), priorsigma = c(1), priorphi = c(5, 1.5), priorbeta =  c(0, 10000))
summary(res[["para"]][[1]][,-c(4,5)])
summary(res[["beta"]])
plot(res)
ht <- res[["latent"]][[1]]
library(fanplot)
require(latex2exp)
df <- as.data.frame(ht)
plot(NULL, main="Percentiles", xlim = c(1, T+1), ylim = c(-12.5, -7.5), xlab = "Time", ylab = TeX("$h_{t}$"))
fan(data = df)
lines(colMeans(ht), col = "black", lw = 2)
lines(h, col = "blue", lw = 2)

########################## Stochastic volatility models: SMC approach ########################## 
rm(list = ls())
set.seed(010101)
T <- 1250; K <- 2
mu <- -10; phi <- 0.95; sigma <- 0.3
h <- numeric(T)
y <- numeric(T)
h[1] <- rnorm(1, mu, sigma / sqrt(1 - phi^2))  # Initial state
y[1] <- rnorm(1, 0, exp(h[1] / 2))           # Initial observation
for (t in 2:T) {
  h[t] <- mu + phi*(h[t-1]-mu) + rnorm(1, 0, sigma)
  y[t] <- rnorm(1, 0, sd = exp(0.5*h[t]))
}

pmhOutput <- pmhtutorial::particleMetropolisHastingsSVmodel(y,
                                               initialTheta = c(0, 0.9, 0.2),
                                               noParticles=500,
                                               noIterations=1000,
                                               stepSize=diag(c(0.05, 0.0002, 0.002)))
summary(coda::mcmc(pmhOutput[["theta"]]))
mean(pmhOutput[["proposedThetaAccepted"]])
ht <- pmhOutput[["xHatFiltered"]]
library(fanplot)
require(latex2exp)
df <- as.data.frame(ht)
plot(NULL, main="Percentiles", xlim = c(1, T+1), ylim = c(-13.5, -6.5), xlab = "Time", ylab = TeX("$h_{t}$"))
fan(data = df)
lines(colMeans(ht), col = "black", lw = 2)
lines(h, col = "blue", lw = 2)

