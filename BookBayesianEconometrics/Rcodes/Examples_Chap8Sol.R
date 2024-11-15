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