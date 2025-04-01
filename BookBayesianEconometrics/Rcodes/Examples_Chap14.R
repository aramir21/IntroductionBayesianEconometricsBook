######### ABC Exchange rate og returns: USD/EURO
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
