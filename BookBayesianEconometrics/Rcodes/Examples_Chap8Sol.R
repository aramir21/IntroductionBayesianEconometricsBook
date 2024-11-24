########################## Simulation exercise: Dynamic linear model ########################## 
rm(list = ls()); set.seed(010101)
T <- 100; sig2 <- 0.5^2; r <- 0.5; sigW2 <- sig2*r
x <- rnorm(T, mean = 1, sd = 0.1*sig2^0.5) 
e <- rnorm(T, mean = 0, sd = sig2^0.5)
w <- rnorm(T, mean = 0, sd = sigW2^0.5)
K <- 1 
Bt <- matrix(NA, T, K); Bt[1] <- 1
yt <- rep(NA, T) 
yt[1] <- x[1]*Bt[1] + e[1]
for(t in 1:T){
  if(t == 1){
    Bt[t,] <- w[t]
  }else{
    Bt[t,] <- Bt[t-1,] + w[t]
  }
  yt[t] <- x[t]%*%Bt[t,] + e[t]
}
# State spece model
ModelReg <- function(par){
  Mod <- dlm::dlmModReg(x, dV = exp(par[1]), dW = exp(par[2]), m0 = rep(1, K),
                        C0 = sigW2*diag(K), addInt = FALSE)
  return(Mod)
}
sig2New <- sig2*25
RegFilter1 <- dlm::dlmFilter(yt, ModelReg(c(sig2New, sig2New*0.5)))
ytfil1 <- RegFilter1[["m"]][-1]
RegFilter2 <- dlm::dlmFilter(yt, ModelReg(c(sig2New, sig2New*0.1)))
ytfil2 <- RegFilter2[["m"]][-1]
Time <- 1:T
df1 <- as.data.frame(cbind(Time,yt, ytfil1, ytfil2))
library(ggplot2)
ggplot(df1, aes(x=Time)) +
  geom_line(aes(y=yt), colour="black", linewidth=1, linetype=1) +
  geom_line(aes(y=ytfil1), colour="green", linewidth=1, alpha=0.9, linetype=2) +
  geom_line(aes(y=ytfil2), color="blue", linewidth=1, alpha=0.9, linetype=3)

########################## Kalman filter recursion ##########################
rm(list = ls()); set.seed(010101)
T <- 100; sig2 <- 1; r <- 0.5; sigW2 <- sig2*r
xt <- rnorm(T, mean = 1, sd = 0.1*sig2^0.5) 
e <- rnorm(T, mean = 0, sd = sig2^0.5)
w <- rnorm(T, mean = 0, sd = sigW2^0.5)
K <- 1 
Bt <- matrix(NA, T, K); Bt[1] <- 0
yt <- rep(NA, T) 
for(t in 1:T){
  if(t == 1){
    Bt[t,] <- w[t]
  }else{
    Bt[t,] <- Bt[t-1,] + w[t]
  }
  yt[t] <- xt[t]%*%Bt[t,] + e[t]
}
# Filtering
KFR <- function(y, x, b1, B1, Omega, sig2){
  e <- as.numeric(y - t(x)%*%b1) # Error
  R <- B1 + Omega
  q <- as.numeric(t(x)%*%R%*%x + sig2)
  K <- (R%*%x)/q
  bt <- b1 + K*e
  Bt <- R - R%*%x%*%t(x)%*%R/q
  Result <- list(bt = bt, Bt = Bt) 
  return(Result)
}
KFresb <- list(); KFresB <- list()
b0 <- 0; B0 <- sigW2
for(t in 1:T){
  KFrest <- KFR(y = yt[t], x = xt[t], b1 = b0, B1 = B0, Omega = sigW2, sig2 = sig2)
  b0 <- KFrest[["bt"]]; B0 <- KFrest[["Bt"]]
  KFresb[[t]] <- b0; KFresB[[t]] <- B0 
}
ModelReg <- function(par){
  Mod <- dlm::dlmModReg(xt, dV = exp(par[1]), dW = exp(par[2]), m0 = rep(0, K),
                        C0 = sigW2*diag(K), addInt = FALSE)
  return(Mod)
}
RegFilter1 <- dlm::dlmFilter(yt, ModelReg(c(sig2, sigW2)))
plot(Bt, type = "l")
lines(unlist(KFresb), type = "l", col = "red")
lines(RegFilter1[["m"]][-1], type = "l", col = "blue")
# Smoothing 
BSR <- function(b, B, slead, Slead, Omega){
  b = bfilter[t]; B = Bfilter[t]; slead = slead; Slead = slead; Omega = sigW2
  R1 <- B + Omega
  s <- b + B%*%solve(R1)%*%(slead - b)
  S <- B - B%*%solve(R1)%*%(R1 - Slead)%*%solve(R1)%*%B
  Result <- list(st = s, St = S) 
  return(Result)
}
bfilter <- unlist(KFresb); Bfilter <- unlist(KFresB)
slead <- bfilter[T]; Slead <- Bfilter[T]
smooth <- rep(slead, T); Smooth <- rep(Slead, T);  
for(t in (T-1):1){
  BSRrest <- BSR(b = bfilter[t], B = Bfilter[t], slead = slead, Slead = Slead, Omega = sigW2)
  slead <- BSRrest[["st"]]; Slead <- BSRrest[["st"]]
  smooth[t] <- slead
}
RegSmoth <- dlm::dlmSmooth(yt, ModelReg(c(sig2, sigW2)))
VarSmooth <- dlm::dlmSvd2var(u = RegSmoth[["U.S"]], RegSmoth[["D.S"]])
SDVarSmoothB2 <- sapply(2:(T+1), function(t){VarSmooth[[t]][K,K]^0.5}) 
plot(Bt, type = "l")
lines(smooth, type = "l", col = "orange")
lines(RegSmoth[["s"]][-1], type = "l", col = "green")

########################## Kalman filter recursion DeJong-Shepard 1995 and Bayesian inference ##########################
rm(list = ls()); set.seed(010101)
T <- 200; sig2 <- 0.1; tau2 <- 5; tau <- tau2^0.5; sigW2 <- sig2/tau2
xt <- rnorm(T, mean = 1, sd = 2*sig2^0.5)
zt <- rnorm(T, mean = 1, sd = 2*sig2^0.5)
alpha <- 2
e <- MASS::mvrnorm (T, mu = c(0, 0), Sigma = sig2*diag(2))
K <- 1
h1 <- matrix(c(1, 0), K + 1, 1); h2 <- matrix(c(0, 1/tau), K + 1, 1); t(h1)%*%h2
e1 <- e%*%h1; e2 <- e%*%h2; var(e1); var(e2)
Bt <- matrix(NA, T, K); Bt[1] <- 0
yt <- rep(NA, T) 
for(t in 1:T){
  if(t == 1){
    Bt[t,] <- e2[t]
  }else{
    Bt[t,] <- Bt[t-1,] + e2[t] 
  }
  yt[t] <- zt[t]*alpha + xt[t]*Bt[t] + e1[t]
}
plot(yt, type = "l")
plot(Bt, type = "l")

KFRDS <- function(y, z, x, b, B, alpha, tau2){
  e <- as.numeric(y - z*alpha - t(x)%*%b) # Error
  q <- as.numeric(t(x)%*%B%*%x + t(h1)%*%h1)
  K <- (B%*%x)/q
  bt <- b + K*e
  h2 <- matrix(c(0, 1/tau2^0.5))
  Bt <- B - B%*%x%*%t(K) + t(h2)%*%h2
  Result <- list(bt = bt, Bt = Bt, et = e, Kt = K, qt = q) 
  return(Result)
}
b0 <- 0; B0 <- sigW2
# y = yt[t]; x = xt[t]; b = b0; B = B0; sig2 = sig2; z = zt[t]; alpha = alpha
KFDSresb <- list(); KFDSresB <- list()
KFDSresE <- list(); KFDSresK <- list(); KFDSresQ <- list() 
for(t in 1:T){
  KFDSrest <- KFRDS(y = yt[t], z = zt[t], x = xt[t], b = b0, B = B0, alpha = alpha, tau2 = tau2)
  b0 <- KFDSrest[["bt"]]; B0 <- KFDSrest[["Bt"]]
  KFDSresb[[t]] <- b0; KFDSresB[[t]] <- B0 
  KFDSresE[[t]] <- KFDSrest[["et"]]; KFDSresK[[t]] <- KFDSrest[["Kt"]]; KFDSresQ[[t]] <- KFDSrest[["qt"]]  
}
ModelReg <- function(par){
  Mod <- dlm::dlmModReg(xt, dV = exp(par[1]), dW = exp(par[2]), m0 = rep(0, K),
                        C0 = sigW2*diag(K), addInt = FALSE)
  return(Mod)
}
RegFilter1 <- dlm::dlmFilter(yt - alpha*zt, ModelReg(c(sig2, sigW2)))

KFR <- function(y, x, b, B, Omega, sig2){
  e <- as.numeric(y - t(x)%*%b) # Error
  R <- B + Omega
  q <- as.numeric(t(x)%*%R%*%x + sig2)
  K <- (R%*%x)/q
  bt <- b + K*e
  Bt <- R - R%*%x%*%t(x)%*%R/q
  Result <- list(bt = bt, Bt = Bt) 
  return(Result)
}
# t <- 1
KFresb <- list(); KFresB <- list()
b0 <- 0; B0 <- sigW2
for(t in 1:T){
  # y = yt[t]; x = xt[t]; b = b0; B = B0; Omega = sigW2; sig2 = sig2
  KFrest <- KFR(y = yt[t] - zt[t]*alpha, x = xt[t], b = b0, B = B0, Omega = sigW2, sig2 = sig2)
  b0 <- KFrest[["bt"]]; B0 <- KFrest[["Bt"]]
  KFresb[[t]] <- b0; KFresB[[t]] <- B0 
}

plot(Bt, type = "l")
lines(unlist(KFDSresb), type = "l", col = "red")
lines(unlist(KFresb), type = "l", col = "green")
lines(RegFilter1[["m"]][-1], type = "l", col = "blue")

BSRDS <- function(r, M, Q, e, K, x, sig2, tau2){
  h2 <- matrix(c(0, 1/tau2^0.5))
  Lambda <- t(h2)%*%h2 
  C <- Lambda - Lambda%*%M%*%t(Lambda)
  xi <- rnorm(1, 0, sd = (sig2*C)^0.5)
  L <- 1 - K*x
  V <- Lambda%*%M%*%L
  rtl <- x*e/Q + t(L)*r - t(V)%*%solve(C)%*%xi
  Mtl <- x%*%t(x)/Q + t(L)%*%M%*%L + t(V)%*%solve(C)%*%V
  eta <- Lambda%*%r + xi
  Result <- list(rtl = rtl, Mtl = Mtl, etat = eta)
  return(Result)
}
et <- unlist(KFDSresE); Qt <- unlist(KFDSresQ); Kt <- unlist(KFDSresK)
rT <- 0; MT <- 0
BSRDSreseta <- rep(0, T); BSRDSresetaMean <- rep(0, T)  
# t <- T - 1
for(t in (T-1):1){
  BSRDSrest <- BSRDS(r = rT, M = MT, Q = Qt[t+1], e = et[t+1], K = Kt[t+1], x = xt[t+1], sig2 = sig2, tau2 = tau2)
  rT <- BSRDSrest[["rtl"]]; MT <- BSRDSrest[["Mtl"]]
  BSRDSreseta[t+1] <- BSRDSrest[["etat"]] 
}
Bt1DS <- rep(0,T)
for(t in 2:T){
  Bt1DS[t] <- Bt1DS[t-1] + BSRDSreseta[t]
}
RegSmoth <- dlm::dlmSmooth(yt - alpha*zt, ModelReg(c(sig2, sigW2)))
plot(Bt, type = "l")
lines(RegSmoth[["s"]][-1], type = "l", col = "purple")
lines(Bt1DS, type = "l", col = "orange")
# Recursion functions
# Filter
KFRDSnew <- function(y, z, x, b, B, alpha, tau2){
  e <- as.numeric(y - z*alpha - t(x)%*%b)
  q <- as.numeric(t(x)%*%B%*%x + t(h1)%*%h1)
  K <- (B%*%x)/q
  bt <- b + K*e
  h2 <- matrix(c(0, 1/tau2^0.5))
  Bt <- B - B%*%x%*%t(K) + t(h2)%*%h2
  Result <- list(bt = bt, Bt = Bt, et = e, Kt = K, qt = q) 
  return(Result)
}
# Smooth
BSRDSnew <- function(r, M, Q, e, K, x, sig2, tau2){
  h2 <- matrix(c(0, 1/tau2^0.5))
  Lambda <- t(h2)%*%h2 
  C <- Lambda - Lambda%*%M%*%t(Lambda)
  xi <- rnorm(1, 0, sd = (sig2*C)^0.5)
  L <- 1 - K*x
  V <- Lambda%*%M%*%L
  rtl <- x*e/Q + t(L)*r - t(V)%*%solve(C)%*%xi
  Mtl <- x%*%t(x)/Q + t(L)%*%M%*%L + t(V)%*%solve(C)%*%V
  eta <- Lambda%*%r + xi
  Result <- list(rtl = rtl, Mtl = Mtl, etat = eta)
  return(Result)
}
# Gibbs functions
PostSig2 <- function(bbt, alpha, tau2){
  an <- T*(K + 1) + alpha0
  term1 <- diff(c(0, bbt))
  term2 <- c(yt - alpha*zt - xt*bbt)
  dn <- delta0 + sum(term1^2)*tau2  + sum(term2^2)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PosAlpha <- function(bbt, sig2){
  An <- solve(A0i + sig2^(-1)*t(zt)%*%zt)
  term <- yt-xt*bbt
  an <- An%*%(A0i%*%a0 + sig2^(-1)*t(zt)%*%term)
  Alpha <- MASS::mvrnorm(1, an, An)
  return(Alpha)
}
PostTau2 <- function(bbt, sig2){
  v1n <- v0 + T
  term1 <- diff(c(0,bbt))
  v2n <- v0 + sig2^(-1)*sum(term1^2)
  tau2 <- rgamma(1, v1n/2, v2n/2)
  return(tau2)
}
# Hyperparameter
alpha0 <- 1; delta0 <- 1; v0 <- 5
a0 <- 0; A0 <- 1; A0i <- 1/A0
# tau2post <- tau2; sig2post <- sig2; alphapost <- alpha
tau2post <- 2; sig2post <- 0.2; alphapost <- 0
# Inference
S <- 2000; burnin <- 500; tot <- S + burnin; thin <- 5
# Posterior draws
sig2s <- rep(NA, tot); tau2s <- rep(NA, tot)
alphas <- rep(NA, tot); betas <- matrix(NA, tot, T)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for (s in 1:tot){
  b0 <- 0; B0 <- sigW2
  KFDSresb <- list(); KFDSresB <- list()
  KFDSresE <- list(); KFDSresK <- list(); KFDSresQ <- list()
  for(t in 1:T){
    KFDSrest <- KFRDSnew(y = yt[t], z = zt[t], x = xt[t], b = b0, B = B0, alpha = alphapost, tau2 = tau2post)
    b0 <- KFDSrest[["bt"]]; B0 <- KFDSrest[["Bt"]]
    KFDSresb[[t]] <- b0; KFDSresB[[t]] <- B0
    KFDSresE[[t]] <- KFDSrest[["et"]]; KFDSresK[[t]] <- KFDSrest[["Kt"]]; KFDSresQ[[t]] <- KFDSrest[["qt"]]
  }
  et <- unlist(KFDSresE); Qt <- unlist(KFDSresQ); Kt <- unlist(KFDSresK)
  rT <- 0; MT <- 0
  BSRDSreseta <- rep(0, T); BSRDSresetaMean <- rep(0, T)
  for(t in (T-1):1){
    BSRDSrest <- BSRDSnew(r = rT, M = MT, Q = Qt[t+1], e = et[t+1], K = Kt[t+1], x = xt[t+1], sig2 = sig2post, tau2 = tau2post)
    rT <- BSRDSrest[["rtl"]]; MT <- BSRDSrest[["Mtl"]]
    BSRDSreseta[t+1] <- BSRDSrest[["etat"]]
  }
  Bt1DSpost <- rep(0,T)
  for(t in 1:T){
    if(t == 1){
      Bt1DSpost[t] <- BSRDSreseta[t]
    }else{
      Bt1DSpost[t] <- Bt1DSpost[t-1] + BSRDSreseta[t]
    }
  }
  # Bt1DSpost <- Bt
  sig2post <- PostSig2(bbt = Bt1DSpost, alpha = alphapost, tau2 = tau2post)
  alphapost <- PosAlpha(bbt = Bt1DSpost, sig2 = sig2post)
  tau2post <- PostTau2(bbt = Bt1DSpost, sig2 = sig2post)
  sig2s[s] <- sig2post
  tau2s[s] <- tau2post
  alphas[s] <- alphapost
  betas[s,] <- Bt1DSpost
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0), "% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
betasF <- betas[keep,]
sig2sF <- coda::mcmc(sig2s[keep])
tau2sF <- coda::mcmc(tau2s[keep])
alphasF <- coda::mcmc(alphas[keep])
summary(sig2sF)
plot(sig2sF)
summary(tau2sF)
plot(tau2sF)
summary(alphasF)
plot(alphasF)
library(fanplot)
df <- as.data.frame(betasF)
plot(NULL, main="Percentiles", xlim = c(1, T+1), ylim = c(-4, 1), xlab = "Time", ylab = TeX("$\\beta_{t1}$"))
fan(data = df)
lines(colMeans(betasF), col = "black", lw = 2)
lines(Bt, col = "blue", lw = 2)

########################## Hamiltonian Monte Carlo: AR(2) model ########################## 
# Simulation AR(2)
rm(list = ls()); set.seed(010101); T <- 1000; K <- 4 
mu <- 0.5; phi1 <- 0.5; phi2 <- 0.3; sig <- 0.5 
Ey <- mu/(1-phi1-phi2); Sigy <- sig*((1-phi2)/(1-phi2-phi1^2-phi2*phi1^2-phi2^2+phi2^3))^0.5 
y <- rnorm(T, mean = Ey, sd = Sigy); e <- rnorm(T, mean = 0, sd = sig)
for(t in 3:T){
  y[t] <- mu + phi1*y[t-1] + phi2*y[t-2] + e[t]
}
# Hyperparameters
d0 <- 0.01; a0 <- 0.01; mu0 <- 0; MU0 <- 1
phi0 <- c(0, 0); Phi0 <- diag(2)
# Log posterior multiply by -1 to use optim
LogPost <- function(theta, y){
  mu <- theta[1]; phi1 <- theta[2]; phi2 <- theta[3]
  tau <- theta[4]; sig2 <- exp(tau); logLik <- NULL
  for(t in 3:T){
    logLikt <- dnorm(y[t], mean = mu + phi1*y[t-1] + phi2*y[t-2], sd = sig2^0.5, log = TRUE)
    logLik <- c(logLik, logLikt)
  }
  logLik <- sum(logLik)
  logPrior <- dnorm(mu, mean = mu0, sd = MU0^0.5, log = TRUE) + dnorm(phi1, mean = phi0[1], sd = Phi0[1,1]^0.5, log = TRUE) + dnorm(phi2, mean = phi0[2], sd = Phi0[2,2]^0.5, log = TRUE) + invgamma::dinvgamma(sig2, shape = a0/2, rate = d0/2, log = TRUE)
  logPosterior <- logLik + logPrior + tau
  return(-logPosterior) # Multiply by -1 to minimize using optim
}
theta0 <- c(mean(y), 0, 0, var(y))
Opt <- optim(theta0, LogPost, y = y, hessian = TRUE)
theta0 <- Opt$par; VarPost <- solve(Opt$hessian)
# Gradient log posterior
GradientTheta <- function(theta, y){
  mu <- theta[1]; phi1 <- theta[2]; phi2 <- theta[3]
  tau <- theta[4]; sig2 <- exp(tau); SumLik <- matrix(0, 3, 1)
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
StatRest <- function(phi1, phi2){
  if(abs(phi2) < 1 & phi1 + phi2 < 1 & phi2 - phi1 < 1){
    check <- 1
  }else{
    check <- 0
  }
  return(check)
}
StatRest(phi1 = -1.5, phi2 = -0.3)
HMC <- function(theta, y, epsilon, M){
  L <- ceiling(1/epsilon)
  Minv <- solve(M)
  thetat <- theta
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
  if(u <= alpha & StatRest(phi1 = theta[2], phi2 = theta[3]) == 1){
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
theta0 <- theta0
M <- solve(VarPost); epsilon0 <- 0.4
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
summary(exp(thetaF[,K]))
ProbAceptF <- coda::mcmc(ProbAcept[keep])
summary(ProbAceptF)

########################## Stochastic volatility models: Sequential importance sampling ########################## 
rm(list = ls())
set.seed(010101)
T <- 250
mu <- -10; phi <- 0.95; sigma <- 0.3
h <- numeric(T)
y <- numeric(T)
h[1] <- rnorm(1, mu, sigma / sqrt(1 - phi^2))  # Initial state
y[1] <- rnorm(1, 0, exp(h[1] / 2))           # Initial observation
for (t in 2:T) {
  h[t] <- mu + phi*(h[t-1]-mu) + rnorm(1, 0, sigma)
  y[t] <- rnorm(1, 0, sd = exp(0.5*h[t]))
}
# Sequential importance sampling with Conditional Prior as Proposal. See SV example and SIS algorithm in https://www.cs.ubc.ca/~arnaud/doucet_johansen_tutorialPF.pdf
N <- 10000
log_Weights <- matrix(NA, N, T)  # Log weights
Weights <- matrix(NA, N, T)  # Weights 
WeightsST <- matrix(NA, N, T)  # Normalized weights 
particles <- matrix(NA, N, T)   # Particles
logalphas <- matrix(NA, N, T)   # Incremental importance weights
ESS <- rep(NA, T) # Effective sample size
# Initialize particles
particles[, 1] <- rnorm(N, mu, sigma / sqrt(1 - phi^2))  # Stationary prior
log_Weights[, 1] <- dnorm(y[1], 0, sd = exp(0.5*particles[,1]), log = TRUE)  # Likelihood
# Normalize weights
# Weights[, 1] <- exp(log_weights[, 1] - max(log_weights[, 1]))
Weights[, 1] <- exp(log_Weights[, 1])
WeightsST[, 1] <- Weights[, 1] / sum(Weights[, 1])
ESS[1] <- (sum(WeightsST[, 1]^2))^(-1)
# Loop over time
for (t in 2:T) {
  # Sample particles
  particles[, t] <- rnorm(N, mu + phi*(particles[, t - 1] - mu), sigma)  # Sample from proposal
  logalphas[, t] <- dnorm(y[t], 0, sd = exp(0.5*particles[,t]), log = TRUE) 
  log_Weights[, t] <- log_Weights[, t - 1] + logalphas[, t]
  Weights[, t] <- exp(log_Weights[, t])
  WeightsST[, t] <- Weights[, t] / sum(Weights[, t])
  ESS[t] <- (sum(WeightsST[, t]^2))^(-1)
}
FilterDist <- colSums(particles * WeightsST)
SDFilterDist <- (colSums(particles^2 * WeightsST) - FilterDist^2)^0.5
MargLik <- colMeans(Weights)
# Plot results
library(dplyr)
library(ggplot2)
require(latex2exp)
ggplot2::theme_set(theme_bw())
df <- tibble(t = seq(1, T),
             mean = FilterDist,
             lower = FilterDist - 1.96*SDFilterDist,
             upper = FilterDist + 1.96*SDFilterDist,
             x_true = h)

plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1,
                fill = "lightblue") +
    geom_line(aes(y = x_true), colour = "black", alpha = 1,
              linewidth = 0.5) +
    geom_line(aes(y = mean), colour = "blue", linewidth = 0.5) +
    ylab(TeX("$h_{t}$")) + xlab("Time")
  print(p)
}

plot_filtering_estimates(df)
plot(ESS, type = "l", ylab = "Effective sample size", xlab = "Time")

########################## Stochastic volatility models: SMC approach ########################## 
rm(list = ls())
set.seed(010101)
T <- 1250; mu <- -10; phi <- 0.95; sigma <- 0.3
h <- numeric(T)
y <- numeric(T)
h[1] <- rnorm(1, mu, sigma / sqrt(1 - phi^2))  # Initial state
y[1] <- rnorm(1, 0, exp(h[1] / 2))           # Initial observation
for (t in 2:T) {
  h[t] <- mu + phi*(h[t-1]-mu) + rnorm(1, 0, sigma)
  y[t] <- rnorm(1, 0, sd = exp(0.5*h[t]))
}
# Sequential Monte Carlo with Conditional Prior as Proposal. See SV example and SMC algorithm in https://www.cs.ubc.ca/~arnaud/doucet_johansen_tutorialPF.pdf
N <- 10000
log_Weights <- matrix(NA, N, T)  # Log weights
Weights <- matrix(NA, N, T)  # Weights 
WeightsST <- matrix(NA, N, T)  # Normalized weights 
WeightsSTT <- matrix(1/N, N, T)  # Normalized weights bar 
particles <- matrix(NA, N, T)   # Particles
particlesT <- matrix(NA, N, T)   # Particles bar
logalphas <- matrix(NA, N, T)   # Incremental importance weights
ESS <- rep(NA, T) # Efective sample size
cutoff <- 0.5
# Initialize particles
particles[, 1] <- rnorm(N, mu, sigma / sqrt(1 - phi^2))  # Stationary prior
log_Weights[, 1] <- dnorm(y[1], 0, sd = exp(0.5*particles[,1]), log = TRUE)  # Likelihood
# Normalize weights
# Weights[, 1] <- exp(log_weights[, 1] - max(log_weights[, 1]))
Weights[, 1] <- exp(log_Weights[, 1])
WeightsST[, 1] <- Weights[, 1] / sum(Weights[, 1])
ESS[1] <- (sum(WeightsST[, 1]^2))^(-1)
ind <- sample(1:N, size = N, replace = TRUE, prob = WeightsST[, 1]) # Resample 
particles[, 1] <- particles[, 1] # Resampled particles
particlesT[, 1] <- particles[ind, 1] # Resampled particles
WeightsST[, 1] <- rep(1/N, N) # Resampled weights
pb <- winProgressBar(title = "progress bar", min = 0, max = T, width = 300)
for (t in 2:T) {
  particles[, t] <- rnorm(N, mu + phi*(particles[, t - 1] - mu), sigma)  # Sample from proposal
  logalphas[, t] <- dnorm(y[t], 0, sd = exp(0.5*particles[,t]), log = TRUE) 
  Weights[, t] <- exp(logalphas[, t])
  WeightsST[, t] <- Weights[, t] / sum(Weights[, t])
  ESS[t] <- (sum(WeightsST[, t]^2))^(-1)
  particlesT[, t] <- particles[, t]
  if(ESS[t] < N*cutoff){
    ind <- sample(1:N, size = N, replace = TRUE, prob = WeightsST[, t])
    particles[, 1:t] <- particles[ind, 1:t]
    WeightsST[, t] <- 1/N
  }else{if(t == T & ESS[t] < N*cutoff)
    ind <- sample(1:N, size = N, replace = TRUE, prob = WeightsST[, t])
  particlesT[, 1:t] <- particles[ind, 1:t]
  }
  setWinProgressBar(pb, t, title=paste( round(t/T*100, 0), "% done"))
}
close(pb)
FilterDist <- colSums(particles * WeightsST)
SDFilterDist <- (colSums(particles^2 * WeightsST) - FilterDist^2)^0.5
FilterDistT <- colSums(particlesT * WeightsSTT)
SDFilterDistT <- (colSums(particlesT^2 * WeightsSTT) - FilterDistT^2)^0.5
MargLik <- colMeans(Weights)
plot(MargLik, type = "l")
plot(ESS, type = "l")
library(dplyr)
library(ggplot2)
require(latex2exp)
ggplot2::theme_set(theme_bw())
Tfig <- 250
keepFig <- 1:Tfig
df <- tibble(t = keepFig,
             mean = FilterDist[keepFig],
             lower = FilterDist[keepFig] - 2*SDFilterDist[keepFig],
             upper = FilterDist[keepFig] + 2*SDFilterDist[keepFig],
             meanT = FilterDistT[keepFig],
             lowerT = FilterDistT[keepFig] - 2*SDFilterDistT[keepFig],
             upperT = FilterDistT[keepFig] + 2*SDFilterDistT[keepFig],
             x_true = h[keepFig])

plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1,
                fill = "lightblue") +
    geom_line(aes(y = x_true), colour = "black", alpha = 1,
              linewidth = 0.5) +
    geom_line(aes(y = mean), colour = "blue", linewidth = 0.5) +
    geom_line(aes(y = meanT), colour = "purple", linewidth = 0.5) +
    ylab(TeX("$h_{t}$")) + xlab("Time")
  print(p)
}
plot_filtering_estimates(df)

########################## Application of stochastic volatility model: Exchange rate ########################## 
rm(list = ls())
set.seed(010101)
DataExcRate <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/17ExcRate.csv", sep = ",", header = TRUE, quote = "")
attach(DataExcRate)
MCMC <- 10000; burnin <- 10000; thin <- 5
y <- USDEUR - mean(USDEUR)
plot(y, type = "l")
res <- stochvol::svsample(y, draws = MCMC, burnin = burnin, thin = thin, priormu = c(0, 100), priorsigma = c(1), priorphi = c(5, 1.5), priorbeta =  c(0, 10000))
summary(res[["para"]][[1]][,-c(4,5)])
plot(res)
ht <- res[["latent"]][[1]]
library(dplyr)
library(ggplot2)
require(latex2exp)
ggplot2::theme_set(theme_bw())
x_means <- colMeans(ht)
x_quantiles <- apply(ht, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
df <- tibble(t = 1:length(y),
             mean = x_means,
             lower = x_quantiles[1, ],
             upper = x_quantiles[2, ])
plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1,
                fill = "lightblue") +
    geom_line(aes(y = mean), colour = "blue", linewidth = 0.5) +
    ylab(TeX("$h_{t}$")) + xlab("Time")
  print(p)
}
plot_filtering_estimates(df)

y <- USDGBP - mean(USDGBP)
plot(y, type = "l")
res <- stochvol::svsample(y, draws = MCMC, burnin = burnin, thin = thin, priormu = c(0, 100), priorsigma = c(1), priorphi = c(5, 1.5), priorbeta =  c(0, 10000))
summary(res[["para"]][[1]][,-c(4,5)])
plot(res)
ht <- res[["latent"]][[1]]
library(dplyr)
library(ggplot2)
require(latex2exp)
ggplot2::theme_set(theme_bw())
x_means <- colMeans(ht)
x_quantiles <- apply(ht, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
df <- tibble(t = 1:length(y),
             mean = x_means,
             lower = x_quantiles[1, ],
             upper = x_quantiles[2, ])
plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1,
                fill = "lightblue") +
    geom_line(aes(y = mean), colour = "blue", linewidth = 0.5) +
    ylab(TeX("$h_{t}$")) + xlab("Time")
  print(p)
}
plot_filtering_estimates(df)

y <- GBPEUR - mean(GBPEUR)
plot(y, type = "l")
res <- stochvol::svsample(y, draws = MCMC, burnin = burnin, thin = thin, priormu = c(0, 100), priorsigma = c(1), priorphi = c(5, 1.5), priorbeta =  c(0, 10000))
summary(res[["para"]][[1]][,-c(4,5)])
plot(res)
ht <- res[["latent"]][[1]]
library(dplyr)
library(ggplot2)
require(latex2exp)
ggplot2::theme_set(theme_bw())
x_means <- colMeans(ht)
x_quantiles <- apply(ht, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
df <- tibble(t = 1:length(y),
             mean = x_means,
             lower = x_quantiles[1, ],
             upper = x_quantiles[2, ])
plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1,
                fill = "lightblue") +
    geom_line(aes(y = mean), colour = "blue", linewidth = 0.5) +
    ylab(TeX("$h_{t}$")) + xlab("Time")
  print(p)
}
plot_filtering_estimates(df)

