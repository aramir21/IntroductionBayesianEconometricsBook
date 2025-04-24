######### ABC Exchange rate log returns: USD/EURO
rm(list = ls()); set.seed(010101)
library(EasyABC)
dfExcRate <- read.csv(file = "https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/ExchangeRate.csv", sep = ",", header = T)
attach(dfExcRate)
str(dfExcRate)
Day <- as.Date(Date, format = "%d/%m/%Y")
n <- length(USDEUR)
# Summary statistics
SumSt <- function(y) {
  Oct <- quantile(y, c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875))
  eta1 <- Oct[6] - Oct[2]
  eta2 <- (Oct[6] + Oct[2] - 2 * Oct[4]) / eta1
  eta3 <- (Oct[7] - Oct[5] + Oct[3] - Oct[1]) / eta1
  autocor <- acf(y, lag = 2, plot = FALSE)
  autocor[["acf"]][2:3]
  Etay <- c(Oct, eta1, eta2, eta3, autocor[["acf"]][2:3])
  return(Etay)
}

# g-and-k distribution
RGKnewSum <- function(par) {
  z <- NULL
  theta <- par[1]; a <- par[2]; b <- par[3]; g <- par[4]; k <- par[5]
  e <- rnorm(n + 1)
  for(t in 2:(n + 1)){
    zt <- e[t] + theta * e[t-1]
    z <- c(z, zt)
  }
  zs <- z / (1 + theta^2)^0.5
  x <- a + b * (1 + 0.8 * (1 - exp(-g * zs)) / (1 + exp(-g * zs))) * (1 + zs^2)^k * zs
  Etaz <- SumSt(x)
  return(Etaz)
}

toy_prior <- list(c("unif",-1,1), c("unif",0,5), c("unif", 0,5), c("unif", -5,5), c("unif", -0.5,5))
sum_stat_obs <- SumSt(USDEUR)
tick <- Sys.time()
ABC_AR <- ABC_rejection(model=RGKnewSum, prior=toy_prior,
                      summary_stat_target = sum_stat_obs, nb_simul=150000, tol = 0.0067,
                      progress_bar = TRUE)
tock <- Sys.time()
tock - tick
PostABCAR <- coda::mcmc(ABC_AR$param)
summary(PostABCAR)
tick <- Sys.time()
ABC_MCMC <- ABC_mcmc(method="Marjoram", model=RGKnewSum,
                   prior=toy_prior, summary_stat_target=sum_stat_obs,
                   n_rec = 100000, progress_bar = TRUE)
tock <- Sys.time()
tock - tick
PostABCMCMC <- coda::mcmc(ABC_MCMC[["param"]][order(ABC_MCMC[["dist"]])[1:1000],])
summary(PostABCMCMC)

tick <- Sys.time()
ABC_SMC <- ABC_sequential(method="Lenormand", model=RGKnewSum,
                        prior=toy_prior, summary_stat_target=sum_stat_obs,
                        nb_simul = 30000, alpha = 0.034, p_acc_min = 0.05,
                        progress_bar = TRUE)
tock <- Sys.time()
tock - tick
PostABCSMC <- coda::mcmc(ABC_SMC[["param"]])
summary(PostABCSMC)
# Figures 
library(ggplot2); library(latex2exp)
Sp <- 1000
df1 <- data.frame(
  Value = c(PostABCAR[1:Sp,1], PostABCMCMC[1:Sp,1], PostABCSMC[1:Sp,1]),
  Distribution = factor(c(rep("AR", Sp), rep("MCMC", Sp), rep("SMC", Sp)))
)

dentheta <- ggplot(df1, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = TeX("Posterior density plot: $theta$"), x = TeX("$theta$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red", "green")) +  theme_minimal() +
  theme(legend.title = element_blank())

df2 <- data.frame(
  Value = c(PostABCAR[1:Sp,4], PostABCMCMC[1:Sp,4], PostABCSMC[1:Sp,4]),
  Distribution = factor(c(rep("AR", Sp), rep("MCMC", Sp), rep("SMC", Sp)))
)

deng <- ggplot(df2, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = "Posterior density plot: g", x = "g", y = "Posterior density") +
  scale_color_manual(values = c("blue", "red", "green")) +  theme_minimal() +
  theme(legend.title = element_blank())

df3 <- data.frame(
  Value = c(PostABCAR[1:Sp,5], PostABCMCMC[1:Sp,5], PostABCSMC[1:Sp,5]),
  Distribution = factor(c(rep("AR", Sp), rep("MCMC", Sp), rep("SMC", Sp)))
)

denk <- ggplot(df3, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = "Posterior density plot: k", x = "k", y = "Posterior density") +
  scale_color_manual(values = c("blue", "red", "green")) +  theme_minimal() +
  theme(legend.title = element_blank())

library(ggpubr)
ggarrange(dentheta, deng, denk, labels = c("A", "B", "C"), ncol = 3, nrow = 1,
          legend = "bottom", common.legend = TRUE)

######## BSL: g-and-k simulation ############
rm(list = ls()); set.seed(010101)
library(BSL)
# Simulate g-and-k data
RGKnew <- function(par) {
  z <- NULL
  theta <- par[1]; a <- par[2]; b <- par[3]; g <- par[4]; k <- par[5]
  e <- rnorm(n + 1)
  for(t in 2:(n + 1)){
    zt <- e[t] + theta * e[t-1]
    z <- c(z, zt)
  }
  zs <- z / (1 + theta^2)^0.5
  x <- a + b * (1 + 0.8 * (1 - exp(-g * zs)) / (1 + exp(-g * zs))) * (1 + zs^2)^k * zs
  return(x)
}
# Summary statistics
SumSt <- function(y) {
  Oct <- quantile(y, c(0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875))
  eta1 <- Oct[6] - Oct[2]
  eta2 <- (Oct[6] + Oct[2] - 2 * Oct[4]) / eta1
  eta3 <- (Oct[7] - Oct[5] + Oct[3] - Oct[1]) / eta1
  autocor <- acf(y, lag = 2, plot = FALSE)
  autocor[["acf"]][2:3]
  Etay <- c(Oct, eta1, eta2, eta3, autocor[["acf"]][2:3])
  return(Etay)
}
# Prior function
LogPrior <- function(par){
  LogPi <- log(par[1] > -1 & par[1] < 1 & par[2] > 0 & par[2] < 5 & par[3] > 0 & par[3] < 5 & par[4] > -5 & par[4] < 5 & par[5] > -0.5 & par[5] < 5)
  return(LogPi)
}
# Population parameters
theta1 <- 0.8; a <- 1; b <- 0.5; g <- -1; k <- 1
parpop <- c(theta1, a, b, g, k)
K <- 5
n <- 500
y <- RGKnew(par = parpop)
# Algorithm parameters
M <- 200 # Number of iterations to calculate mu and sigma
S <- 11000 # Number of MCMC iterations
burnin <- 1000 # Burn in iterations
thin <- 10 # Thining parameter
keep <- seq(burnin + 1, S, thin)
par0 <- c(0.5, 2, 1, 0, 1) 
Modelgk <- newModel(fnSim = RGKnew, fnSum = SumSt, theta0 = par0, fnLogPrior = LogPrior, verbose = FALSE)
validObject(Modelgk)
# Check if the summary statistics are roughly normal
simgk <- simulation(Modelgk, n = M, theta = par0, seed = 10)
par(mfrow = c(4, 3))
for (i in 1:12){
  eval <- seq(min(simgk$ssx[, i]), max(simgk$ssx[, i]), 0.001)
  densnorm <- dnorm(eval, mean = mean(simgk$ssx[, i]), sd(simgk$ssx[, i])) 
  plot(density(simgk$ssx[, i]), main = "", xlab = "")
  lines(eval, densnorm, col = "red")
}
##### Program from scratch #######
Lims <- matrix(c(-1, 0, 0, -5, -0.5, 1, rep(5, 4)), 5, 2)
parPost <- matrix(NA, S, K)
parPost[1,] <- par0 # arpop # par0
EtaY <- SumSt(y = y)
Zsl <- replicate(20, RGKnew(par = parPost[1,]))
EtasZsl <- t(apply(Zsl, 2, SumSt))
Usl <- colMeans(EtasZsl)
SIGMAsl <- var(EtasZsl)
dnormsl <- mvtnorm::dmvnorm(EtaY, Usl, SIGMAsl, log = TRUE)
tune <- 0.005 # Tuning parameter RW MH
CoVarRW <- tune*diag(K)
accept <- rep(0, S)
tick1 <- Sys.time()
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
for (s in 2:S){
  parc <- MASS::mvrnorm(1, mu = parPost[s-1,], Sigma = CoVarRW)
  RestCheck <- NULL
  for(j in 1:K){
    if(parc[j] < Lims[j,1] | parc[j] > Lims[j,2]){
      Rej <- 1
    }else{
      Rej <- 0
    }
    RestCheck <- c(RestCheck, Rej)
  }
  if(sum(RestCheck) != 0){
    parPost[s,] <- parPost[s-1,]
    accept[s] <- 0
  }else{
    Z <- replicate(M, RGKnew(par = parc))
    EtasZ <- t(apply(Z, 2, SumSt))
    Um <- colMeans(EtasZ)
    SIGMAm <- var(EtasZ)
    dnormc <- mvtnorm::dmvnorm(EtaY, Um, SIGMAm, log = TRUE)
    if(s > round(0.1*S, 0) & s < round(0.2*S, 0)){
      CoVarRW <- var(parPost, na.rm = TRUE)
    }
    if(dnormc == -Inf){
      parPost[s,] <- parPost[s-1,]
      accept[s] <- 0
    }else{
      alpha <- min(1, exp(dnormc - dnormsl))
      U <- runif(1)
      if(U <= alpha){
        parPost[s,] <- parc
        Usl <- Um
        SIGMAsl <- SIGMAm
        dnormsl <- dnormc
        accept[s] <- 1
      }else{
        parPost[s,] <- parPost[s-1,]
        accept[s] <- 0
      }
    }
  }
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
tock1 <- Sys.time()
tock1 - tick1
mean(accept)
PostChainOwn <- coda::mcmc(parPost[keep,])
summary(PostChainOwn)
plot(PostChainOwn)
#### BSL package
tick <- Sys.time()
Resultsgk <- bsl(y = y, n = M, M = S, model = Modelgk, covRandWalk = CoVarRW,
                  method = "BSL", thetaNames = expression(theta, a, b, g, k), plotOnTheFly = TRUE)
tock <- Sys.time()
tock - tick
PostChain <- coda::mcmc(Resultsgk@theta[keep,])
summary(PostChain)
plot(PostChain)
Resultsgk@acceptanceRate
plot(Resultsgk@loglike[keep], type = "l")
sd(Resultsgk@loglike[keep])
# Figures 
library(ggplot2); library(latex2exp)
Sp <- length(keep)
df1 <- data.frame(Value = c(PostChain[,1], PostChainOwn[,1]),
  Distribution = factor(c(rep("BSL", Sp), rep("BSLscratch", Sp))))
dentheta <- ggplot(df1, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = TeX("Posterior density plot: $theta$"), x = TeX("$theta$"), y = "Posterior density") +
  geom_vline(xintercept = theta1, linetype = "dashed", color = "red", linewidth = 1) +
  scale_color_manual(values = c("blue", "red")) +  theme_minimal() +
  theme(legend.title = element_blank())
df2 <- data.frame(Value = c(PostChain[,4], PostChainOwn[,4]),
  Distribution = factor(c(rep("BSL", Sp), rep("BSLscratch", Sp))))
deng <- ggplot(df2, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = "Posterior density plot: g", x = "g", y = "Posterior density") +
  geom_vline(xintercept = g, linetype = "dashed", color = "red", linewidth = 1) +
  scale_color_manual(values = c("blue", "red")) +  theme_minimal() +
  theme(legend.title = element_blank())
df3 <- data.frame(Value = c(PostChain[,5], PostChainOwn[,5]),
  Distribution = factor(c(rep("BSL", Sp), rep("BSLscratch", Sp))))
denk <- ggplot(df3, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = "Posterior density plot: k", x = "k", y = "Posterior density") +
  geom_vline(xintercept = k, linetype = "dashed", color = "red", linewidth = 1) +
  scale_color_manual(values = c("blue", "red")) +  theme_minimal() +
  theme(legend.title = element_blank())

library(ggpubr)
ggarrange(dentheta, deng, denk, labels = c("A", "B", "C"), ncol = 3, nrow = 1,
          legend = "bottom", common.legend = TRUE)

######################### INLA ##########################
# Install Rtools according to your R version: https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html
# Check Rtolls is properly installed using Sys.which("make")
# You should get something like:
# make 
# "C:\\rtools45\\usr\\bin\\make.exe" 
# Install INLA
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
rm(list = ls())
library(INLA)
set.seed(010101)
n <- 10000
x <- rnorm(n, sd = 1)
u <- rnorm(n, sd = 0.5)
intercept <- 1; beta <- 1
id <- 1:n
y <- rpois(n, lambda = exp(intercept + beta * x + u))
my.data <- data.frame(y, x, id)
formula <- y ~ 1 + x + f(id, model="iid")
inla.sla <- inla(formula, data = my.data, family = "poisson", control.compute=list(return.marginals.predictor=TRUE))
inla.ga <- inla(formula, data = my.data, family = "poisson",
                control.inla = list(strategy = "gaussian", int.strategy = "eb"),
                control.compute=list(return.marginals.predictor=TRUE))
inla.la <- inla(formula, data = my.data, family = "poisson",
                control.inla = list(strategy = "laplace",
                                    int.strategy = "grid", dz=0.1, diff.logdens=20),
                control.compute=list(return.marginals.predictor=TRUE))
summary(inla.sla)
marg_sla <- inla.sla$marginals.fixed$x
marg_ga  <- inla.ga$marginals.fixed$x
marg_la  <- inla.la$marginals.fixed$x
plot(marg_sla, type = "l", col = "blue", lwd = 2,
     xlab = expression(beta[x]), ylab = "Density",
     main = "Posterior of βₓ under different INLA strategies")
lines(marg_ga, col = "green", lwd = 2, lty = 2)
lines(marg_la, col = "red", lwd = 2, lty = 3)
abline(v = 1, col = "black", lty = 4, lwd = 2)
legend("topright",
       legend = c("Simplified Laplace", "Gaussian", "Full Laplace", "True βₓ = 1"),
       col = c("blue", "green", "red", "black"),
       lwd = 2, lty = c(1, 2, 3, 4), bty = "n")
# Summary beta sla
inla.qmarginal(c(0.025, 0.5, 0.975), marg_sla) # 95% credible interval

# Variance
marg.prec.sla <- inla.sla$marginals.hyperpar[["Precision for id"]]
marg.prec.ga  <- inla.ga$marginals.hyperpar[["Precision for id"]]
marg.prec.la  <- inla.la$marginals.hyperpar[["Precision for id"]]
marg.var.sla <- inla.tmarginal(function(x) 1 / x, marg.prec.sla)
marg.var.ga  <- inla.tmarginal(function(x) 1 / x, marg.prec.ga)
marg.var.la  <- inla.tmarginal(function(x) 1 / x, marg.prec.la)
# Base plot
plot(marg.var.sla, type = "l", col = "blue", lwd = 2,
     xlab = expression(sigma^2), ylab = "Density",
     main = "Posterior of Random Effect Variance")
lines(marg.var.ga, col = "green", lwd = 2, lty = 2)
lines(marg.var.la, col = "red", lwd = 2, lty = 3)
abline(v = 0.25, col = "black", lty = 4, lwd = 2)
legend("topright", legend = c("Simplified Laplace", "Gaussian", "Full Laplace", "True Variance (0.25)"),
       col = c("blue", "green", "red", "black"), lwd = 2, lty = c(1, 2, 3, 4), bty = "n")
# Summary variance sla
inla.qmarginal(c(0.025, 0.5, 0.975), marg.var.sla) # 95% credible interval

######### Variational Bayes: Linear regression ##########
set.seed(010101)
library(LaplacesDemon)
N <- 500; K <- 2
sig2 <- 1; B <- rep(1, K + 1)
X <- cbind(1, matrix(rnorm(N*K), N, K))
e <- rnorm(N, 0, sig2^0.5)
y <- X%*%B + e
######################### Data List Preparation #########################
mon.names <- "mu[1]"
parm.names <- as.parm.names(list(beta=rep(0,K + 1), sigma2=0))
pos.beta <- grep("beta", parm.names)
pos.sigma2 <- grep("sigma2", parm.names)
PGF <- function(Data) {
  beta <- rnorm(Data$K)
  sigma2 <- runif(1)
  return(c(beta, sigma2))
}
MyData <- list(K=K, PGF=PGF, X=X, mon.names=mon.names,
               parm.names=parm.names, pos.beta=pos.beta, pos.sigma2=pos.sigma2, y=y)
########################## Model Specification ##########################
b0 <- 0; B0 <- 1000; a0 <- 0.01; d0 <- 0.01
Model <- function(parm, Data)
{
  ### Parameters
  beta <- parm[Data$pos.beta]
  sigma2 <- interval(parm[Data$pos.sigma2], 1e-100, Inf)
  parm[Data$pos.sigma2] <- sigma2
  ### Log-Prior
  beta.prior <- sum(dnormv(beta, b0, B0, log=TRUE))
  sigma2.prior <- dinvgamma(sigma2, a0/2, d0/2, log=TRUE)
  ### Log-Likelihood
  mu <- tcrossprod(Data$X, t(beta))
  LL <- sum(dnorm(Data$y, mu, sigma2^0.5, log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + sigma2.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=mu[1],
                   yhat=rnorm(length(mu), mu, sigma2^0.5), parm=parm)
  return(Modelout)
}
Initial.Values <- rep(0,K+2); S <- 10000
Fit <- VariationalBayes(Model, Initial.Values, Data=MyData, Covar=NULL, 
                        Iterations=S, Method="Salimans2", Stop.Tolerance=1e-2, CPUs=1)
print(Fit)
PosteriorChecks(Fit)
caterpillar.plot(Fit, Parms="beta")
plot(Fit, MyData, PDF=FALSE)
Pred <- predict(Fit, Model, MyData, CPUs=1)
summary(Pred, Discrep="Chi-Square")
####### MCMC #######
# Posterior parameters
b0 <- rep(b0, K+1); B0 <- B0*diag(K+1)
bhat <- solve(t(X)%*%X)%*%t(X)%*%y
Bn <- as.matrix(Matrix::forceSymmetric(solve(solve(B0) + t(X)%*%X))) # Force this matrix to be symmetric
bn <- Bn%*%(solve(B0)%*%b0 + t(X)%*%X%*%bhat)
dn <- as.numeric(d0 + t(y)%*%y+t(b0)%*%solve(B0)%*%b0-t(bn)%*%solve(Bn)%*%bn)
an <- a0 + N
# Posterior draws
sig2 <- MCMCpack::rinvgamma(S,an/2,dn/2)
summary(coda::mcmc(sig2))
Betas <- t(sapply(1:S, function(s){MASS::mvrnorm(1, bn, sig2[s]*Bn)}))
summary(coda::mcmc(Betas))

####### VB from scratch #######
dnVB <- ((a0+N+K)/(a0+N))*dn 
anVB <- a0 + N + K
BnVB <- Bn
bnVB <- bn
sig2VB <- MCMCpack::rinvgamma(S,anVB/2,dnVB/2) 
BetasVB <- MASS::mvrnorm(S, mu = bnVB, Sigma = (dn/an)*BnVB)
summary(coda::mcmc(sig2VB))
summary(coda::mcmc(BetasVB))
ELBO <- -N/2*log(2*pi) + a0/2*log(d0/2) - anVB/2*log(dnVB/2) - 0.5*log(det(B0)) + 
  0.5*log(det(BnVB)) - lgamma(a0/2) + lgamma(anVB/2) - K/2*log(anVB/dnVB) + K/2

LogMarLik <- -N/2*log(2*pi) + a0/2*log(d0/2) - an/2*log(dn/2) - 0.5*log(det(B0)) + 
  0.5*log(det(Bn)) - lgamma(a0/2) + lgamma(an/2) 
ELBO; LogMarLik; ELBO < LogMarLik
# CAVI
ELBOfunc <- function(d,B){
  ELBOi <- -N/2*log(2*pi) + a0/2*log(d0/2) - anVB/2*log(d/2) - 0.5*log(det(B0)) + 
    0.5*log(det(B)) - lgamma(a0/2) + lgamma(anVB/2) - K/2*log(anVB/d) + 0.5*(anVB/d)*tr(B%*%solve(Bn))
  return(ELBOi)
}
d <- 100; B <- diag(K) 
ELBOfunc(d = d, B = Bn)
Esig2inv <- anVB/d
ELBOs <- rep(-Inf, S)
epsilon <- 1e-5
for(s in 2:S){
  B <- Esig2inv^(-1)*Bn
  EbQb <- sum(diag(B%*%solve(Bn))) - t(bn)%*%solve(Bn)%*%bn + t(y)%*%y + t(b0)%*%solve(B0)%*%b0
  d <- EbQb + d0
  Esig2inv <- as.numeric(anVB/d)
  ELBOs[s] <- ELBOfunc(d = d, B = B)
  # Check if lower bound decreases
  if (ELBOs[s] < ELBOs[s - 1]) { message("Lower bound decreases!\n")}
  # Check for convergence
  if (ELBOs[s] - ELBOs[s - 1] < epsilon) { break }
  # Check if VB converged in the given maximum iterations
  if (s == S) {warning("VB did not converge!\n")}
}
sig2VBscratch <- MCMCpack::rinvgamma(S,anVB/2,d/2) 
BetasVBscratch <- MASS::mvrnorm(S, mu = bnVB, Sigma = B)
summary(coda::mcmc(sig2VBscratch))
summary(coda::mcmc(BetasVBscratch))

