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

########################## Kalman filter recursion DeJong-Shepard 1995 ##########################
rm(list = ls()); set.seed(010101)
T <- 100; sig2 <- 4; r <- 0.25; sigW2 <- sig2*r
xt <- rnorm(T, mean = 1, sd = 2*sig2^0.5)
zt <- rnorm(T, mean = 1, sd = 1*sig2^0.5)
alpha <- 0.5
e <- MASS::mvrnorm (T, mu = c(0, 0), Sigma = sig2*diag(2))
K <- 1
h1 <- matrix(c(1, 0), K + 1, 1); h2 <- matrix(c(0, r), K + 1, 1); t(h1)%*%h2
e1 <- e%*%h1; e2 <- e%*%h2
Bt <- matrix(NA, T+1, K); Bt[1] <- 0
yt <- rep(NA, T) 
yt[1] <- zt[1]*alpha + xt[1]*Bt[1] + e1[1]
for(t in 2:T){
  Bt[t,] <- Bt[t-1,] + e2[t]
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
t <- 1
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








library(fanplot)
df <- as.data.frame(B2t)
plot(NULL, main="Percentiles", xlim = c(1, T+1), ylim = c(-4, 8), xlab = "Time", ylab = TeX("$\\beta_{t1}$"))
fan(data = df)
lines(Bt[,2], col = "black", lw = 2) 


library(fanplot)
# empty plot
plot(NULL, main="Percentiles", xlim = c(1, T), ylim = c(-1, 7), xlab = "Time", ylab = )
df <- as.data.frame(cbind(Time, SmoothB2, LimInfB2, LimSupB2))
# add fan
fan(data = df)

plot(net, xlim=c(1975,2020), ylim=c(-100,300))
fan(mm, type="interval", start=2013, 
    anchor=net[time(net)==2012])




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
