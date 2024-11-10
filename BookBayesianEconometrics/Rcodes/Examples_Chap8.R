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
gibbsOut <- dlm::dlmGibbsDIG(yt, mod = dlm::dlmModReg(Xt),
                            shape.y = 0.1, rate.y = 0.1,
                            shape.theta = 0.1, rate.theta = 0.1,
                            n.sample = MCMC,
                            thin = thin, save.states = TRUE)

B2t <- matrix(0, MCMC - burnin, T + 1)
for(t in 1:(T+1)){
  B2t[,t] <- gibbsOut[["theta"]][t,3,-c(1:burnin)] 
}

Lims <- apply(B2t, 2, function(x){quantile(x, c(0.025, 0.975))})
# Figure
require(latex2exp) # LaTeX equations in figures
xx <- c(1:(T+1), (T+1):1)
yy <- c(Lims[1,], rev(Lims[2,]))
plot   (xx, yy, type = "n", xlab = "Time", ylab = TeX("$\\beta_{t1}$"))
polygon(xx, yy, col = "lightblue", border = "lightblue")
lines(colMeans(B2t), col = "red", lw = 2)
title("State vector: Inflation rate coefficient")

dV <- coda::mcmc(gibbsOut[["dV"]][-c(1:burnin)])
dW <- coda::mcmc(gibbsOut[["dW"]][-c(1:burnin),])
summary(dV)
summary(dW)
plot(dV)
plot(dW)

library(fanplot)
df <- as.data.frame(B2t)
plot(NULL, main="Percentiles", xlim = c(1, T+1), ylim = c(-2, 1), xlab = "Time", ylab = TeX("$\\beta_{t1}$"))
fan(data = df)
lines(colMeans(B2t), col = "black", lw = 2)
abline(h=0, col = "blue")








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
