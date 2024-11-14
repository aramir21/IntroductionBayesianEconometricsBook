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

########################## Application: Dynamic linear model ########################## 
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
df <- as.data.frame(B2t)
plot(NULL, main="Percentiles", xlim = c(1, T+1), ylim = c(-1, 2), xlab = "Time", ylab = TeX("$\\beta_{t2}$"))
fan(data = df)
lines(colMeans(B2t), col = "black", lw = 2)
abline(h=0, col = "blue")

########################## Kalman filter recursion DeJong-Shepard 1995 and Bayesian inference ##########################
rm(list = ls()); set.seed(010101)
T <- 200; sig2 <- 1; tau <- 2; sigW2 <- sig2/tau^2
xt <- rnorm(T, mean = 1, sd = 2*sig2^0.5)
zt <- rnorm(T, mean = 1, sd = 1*sig2^0.5)
alpha <- 0.5
e <- MASS::mvrnorm (T, mu = c(0, 0), Sigma = sig2*diag(2))
K <- 1
h1 <- matrix(c(1, 0), K + 1, 1); h2 <- matrix(c(0, 1/tau), K + 1, 1); t(h1)%*%h2
e1 <- e%*%h1; e2 <- e%*%h2; var(e1); var(e2)
Bt <- matrix(NA, T, K); Bt[1] <- 0
yt <- rep(NA, T) 
yt[1] <- zt[1]*alpha + xt[1]*Bt[1] + e1[1]
for(t in 1:T){
  if(t == 1){
    Bt[t,] <- e2[t]
  }else{
    Bt[t,] <- Bt[t-1,] + e2[t] 
  }
  yt[t] <- zt[t]*alpha + xt[t]*Bt[t] + e1[t]
}
plot(Bt, type = "l")
plot(yt, type = "l")

KFRDS <- function(y, z, x, b, B, alpha){
  e <- as.numeric(y - z*alpha - t(x)%*%b) # Error
  q <- as.numeric(t(x)%*%B%*%x + t(h1)%*%h1)
  K <- (B%*%x)/q
  bt <- b + K*e
  Bt <- B - B%*%x%*%t(K) + t(h2)%*%h2
  Result <- list(bt = bt, Bt = Bt, et = e, Kt = K, qt = q) 
  return(Result)
}
b0 <- 0; B0 <- sigW2
# y = yt[t]; x = xt[t]; b = b0; B = B0; sig2 = sig2; z = zt[t]; alpha = alpha
KFDSresb <- list(); KFDSresB <- list()
KFDSresE <- list(); KFDSresK <- list(); KFDSresQ <- list() 
for(t in 1:T){
  KFDSrest <- KFRDS(y = yt[t], z = zt[t], x = xt[t], b = b0, B = B0, alpha = alpha)
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

BSRDS <- function(r, M, Q, e, K, x){
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
  BSRDSrest <- BSRDS(r = rT, M = MT, Q = Qt[t+1], e = et[t+1], K = Kt[t+1], x = xt[t+1])
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
PostSig2 <- function(bbt, alpha, tau2){
  an <- T*(K + 1) + alpha0
  term1 <- diff(bbt)
  term2 <- c(yt - alpha*zt - xt*bbt)
  dn <- delta0 + sum(term1^2)*tau2  + sum(term2^2)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostSig2(bbt = Bt, alpha = alpha, tau2 = tau^2)
PosAlpha <- function(bbt, sig2, tau2){
  An <- solve(A0i + sig2^(-1)*t(zt)%*%zt)
  term <- yt-xt*bbt
  an <- An%*%(A0i%*%a0 + sig2^(-1)*t(zt)%*%term)
  Alpha <- MASS::mvrnorm(1, an, An)
  return(Alpha)
}
# Hyperparameter
alpha0 <- 1; delta0 <- 1
a0 <- 0; A0 <- 1; A0i <- 1/A0
PosAlpha(bbt = Bt, sig2 = sig2, tau2 = tau^2)






par01 <- log(SumRegLS$sigma^2); par023 <- log(diag(VarBp))
outMLEReg <- dlm::dlmMLE(yt, parm = c(par01, par023), ModelReg)
exp(outMLEReg$par)
RegSmoth <- dlm::dlmSmooth(yt, ModelReg(outMLEReg$par))
SmoothB2 <- RegSmoth$s[-1,2]
VarSmooth <- dlm::dlmSvd2var(u = RegSmoth[["U.S"]], RegSmoth[["D.S"]])
SDVarSmoothB2 <- sapply(2:(T+1), function(t){VarSmooth[[t]][2,2]^0.5}) 
LimInfB2 <- SmoothB2 - qnorm(0.975)*SDVarSmoothB2
LimSupB2 <- SmoothB2 + qnorm(0.975)*SDVarSmoothB2
# Figure
require(latex2exp) # LaTeX equations in figures
xx <- c(1:T, T:1)
yy <- c(LimInfB2, rev(LimSupB2))
plot   (xx, yy, type = "n", xlab = "Time", ylab = TeX("$\\beta_{t1}$"))
polygon(xx, yy, col = "blue", border = "red")
lines(SmoothB2, col = "red", lw = 2)
title("State vector: Slope parameter")


library(fanplot)
# empty plot
plot(NULL, main="Percentiles", xlim = c(1, T), ylim = c(-3, 7))
df <- as.data.frame(cbind(Time, SmoothB2, LimInfB2, LimSupB2))
# add fan
fan(data = df)




; phi1 <- 0.5
K <- length(B)
e <- rnorm(T, mean = 0, sd = sig2^0.5)
mu <- rnorm(T, mean = 0, sd = (sig2/(1-phi1^2))^0.5)
for(t in 2:T){
  mu[t] <- phi1*mu[t-1] + e[t]
}
y <- X%*%B + mu
plot(y, type = "l")
OLSreg <- lm(y  ~ x)
SumOLSreg <- summary(OLSreg)
SumOLSreg; SumOLSreg$sigma
acf(OLSreg$residuals) 
pacf(OLSreg$residuals)
SSreg <- dlm::dlmModReg(x, dV = SumOLSreg$sigma^2, dW = rep(0, K), m0 = rep(0, K),
                   C0 = diag(K))
RegFilter <- dlm::dlmFilter(y, SSreg)
RegFilter[["m"]][T+1,] # Posterior means coefficients
# State spece model
ModelReg <- function(par){
  Mod <- dlm::dlmModReg(x, dV = exp(par[1]), dW = exp(par[2:3]), m0 = rep(0, K),
                 C0 = diag(K))
  return(Mod)
}
outMLEReg <- dlm::dlmMLE(y, parm = rep(0, 3), ModelReg)
exp(outMLEReg$par)
RegSmoth <- dlm::dlmSmooth(y, ModelReg(outMLEReg$par))
plot(ts(RegSmoth$s[-1,1]))
plot(ts(RegSmoth$s[-1,2]))


########################## Simulation exercise: AR(2) model ########################## 
rm(list = ls())
set.seed(010101)
T <- 150
mu <- 0 # 0.007; 
phi1 <- 0.5; phi2 <- 0.3; sig <- 0.035
Ey <- mu/(1-phi1-phi2); Sigy <- sig*((1-phi2)/(1-phi2-phi1^2-phi2*phi1^2-phi2^2+phi2^3))^0.5 
y <- rnorm(T, mean = Ey, sd = Sigy)
e <- rnorm(T, mean = 0, sd = sig)
for(t in 3:T){
  y[t] <- mu + phi1*y[t-1] + phi2*y[t-2] + e[t]
}
mean(y); sd(y)
y <- ts(y, start=c(1920, 1), frequency=1)
plot(y)
iter <- 1000; burnin <- 500
library(bayesforecast)
sf1 <- bayesforecast::stan_sarima(y,order = c(2, 0, 0), prior_mu0 = normal(0, 1),
                                  prior_ar = normal(0, 1), prior_sigma0 = inverse.gamma(0.01, 0.01),
                  seasonal = c(0, 0, 0), iter = iter, warmup = burnin, chains = 1)
Postmu <- sf1[["stanfit"]]@sim[["samples"]][[1]][["mu0"]][-c(1:burnin)]
Postsig <- sf1[["stanfit"]]@sim[["samples"]][[1]][["sigma0"]][-c(1:burnin)]
Postphi1 <- sf1[["stanfit"]]@sim[["samples"]][[1]][["ar0[1]"]][-c(1:burnin)]
Postphi2 <- sf1[["stanfit"]]@sim[["samples"]][[1]][["ar0[2]"]][-c(1:burnin)]
Postdraws <- cbind(Postmu, Postsig, Postphi1, Postphi2)
summary(coda::mcmc(Postdraws))
# print(sf1)


modelAR2 <- function(par){
  AR2 <- dlm::dlmModARMA(ar = par[1:2], sigma2 = exp(par[3]))
  return(AR2)
}
init <- c(0, 0, 1.1)
parMLE <- dlm::dlmMLE(y, init, modelAR2)
parMLE$par; (exp(parMLE$par[3]))^0.5
dlmAR2 <- modelAR2(parMLE$par)
AR2FilterAR2 <- dlm::dlmFilter(y, dlmAR2)
AR2SmoothAR2 <- dlm::dlmSmooth(AR2FilterAR2)
ts.plot(cbind(y, AR2SmoothAR2$s[-1, 1]), col = 1 : 2)
cbind(y, AR2FilterAR2[["m"]])
