############ g-and-k distribution using ABC ############
rm(list = ls()); set.seed(010101)
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
# Population parameters
theta1 <- 0.8; a <- 1; b <- 0.5; g <- -1; k <- 1
parpop <- c(theta1, a, b, g, k)
n <- 500
y <- RGKnew(par = parpop) 
plot(y, type = "l")
##### ABC Function#####
ABC <- function(S, a, y) {
  prior <- cbind(runif(S,-1,1), runif(S,0,5), runif(S,0,5), runif(S,-5,5), runif(S,-0.5,5))
  Z <- apply(prior, 1, RGKnew)
  EtasZ <- apply(Z, 2, SumSt)
  Etay <- SumSt(y) 
  Dist <- sapply(1:S, function(l) {
    dist(rbind(Etay, EtasZ[, l]))
  })
  OrdPrior <- prior[order(Dist), ]
  SelPrior <- OrdPrior[1:round(S * a), ]
  SelSumSt <- t(EtasZ)[1:round(S * a), ]
  return(list(SelPrior = SelPrior, SelSumSt = SelSumSt))
}
S <- 1000000
a <- 0.001
tick <- Sys.time()
ResABC <- ABC(S = S, a = 0.001, y = y)
tock <- Sys.time()
tock - tick
PostABC_ARown <- ResABC[["SelPrior"]]
# Regression adjusted ABC
X <- ResABC[["SelSumSt"]]-matrix(SumSt(y), S*a, 12, byrow = TRUE)
PostABC_ARownRegAd <- PostABC_ARown
for(j in 1:5){
  Reg <- lm(PostABC_ARown[,j] ~ X)
  # Coefficient of regressor 9 is na.
  PostABC_ARownRegAd[,j] <- PostABC_ARown[,j] - X[,-9]%*%Reg$coefficients[-c(1,9)]
}

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
sum_stat_obs <- SumSt(y)
toy_prior <- list(c("unif",-1,1), c("unif",0,5), c("unif", 0,5), c("unif", -5,5), c("unif", -0.5,5))
library(EasyABC)
tick <- Sys.time()
ABC_AR <- ABC_rejection(model=RGKnewSum, prior=toy_prior,
                        summary_stat_target = sum_stat_obs, nb_simul=260000, tol = 0.00385,
                        progress_bar = TRUE)
tock <- Sys.time()
tock - tick
PostABC_AR <- coda::mcmc(ABC_AR$param)

library(ggplot2); library(latex2exp)
Sp <- 1000
df1 <- data.frame(
  Value = c(PostABC_AR[1:Sp,1], PostABC_ARown[1:Sp,1], PostABC_ARownRegAd[1:Sp,1]),
  Distribution = factor(c(rep("EasyABC", Sp), rep("ABC", Sp), rep("ABCAdj", Sp))))

dentheta <- ggplot(df1, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) + 
  geom_vline(xintercept = theta1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = TeX("Posterior density plot: $theta$"), x = TeX("$theta$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red", "green")) +  theme_minimal() +
  theme(legend.title = element_blank())

df2 <- data.frame(
  Value = c(PostABC_AR[1:Sp,4], PostABC_ARown[1:Sp,4], PostABC_ARownRegAd[1:Sp,4]),
  Distribution = factor(c(rep("EasyABC", Sp), rep("ABC", Sp), rep("ABCAdj", Sp))))

deng <- ggplot(df2, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) + 
  geom_vline(xintercept = g, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = TeX("Posterior density plot: g"), x = TeX("$g$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red", "green")) +  theme_minimal() +
  theme(legend.title = element_blank())

df3 <- data.frame(
  Value = c(PostABC_AR[1:Sp,5], PostABC_ARown[1:Sp,5], PostABC_ARownRegAd[1:Sp,5]),
  Distribution = factor(c(rep("EasyABC", Sp), rep("ABC", Sp), rep("ABCAdj", Sp))))

denk <- ggplot(df3, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +
  geom_vline(xintercept = k, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = TeX("Posterior density plot: k"), x = TeX("$k$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red", "green")) +  theme_minimal() +
  theme(legend.title = element_blank())

library(ggpubr)
ggarrange(dentheta, deng, denk, labels = c("A", "B", "C"), ncol = 3, nrow = 1,
          legend = "bottom", common.legend = TRUE)

summary(coda::mcmc(PostABC_ARown))
summary(coda::mcmc(PostABC_ARownRegAd))
summary(coda::mcmc(PostABC_AR))

#### Adding results from abc package. It gives similar results to EasyABC
library(abc)
S <- 100000
a <- 0.01
prior <- cbind(runif(S,-1,1), runif(S,0,5), runif(S,0,5), runif(S,-5,5), runif(S,-0.5,5))
Z <- apply(prior, 1, RGKnew)
EtasZ <- apply(Z, 2, SumSt)
ABC_ARabc <- abc(target = sum_stat_obs, param = prior, sumstat = t(EtasZ), tol = a,
                 method = "loclinear")
PostABC_ARabc <- ABC_ARabc[["unadj.values"]]
PostABC_ARabc <- ABC_ARabc[["adj.values"]]
# Figures 
library(ggplot2); library(latex2exp)
Sp <- 1000
df1 <- data.frame(
  Value = c(PostABC_AR[1:Sp,1], PostABC_ARown[1:Sp,1], PostABC_ARownRegAd[1:Sp,1], PostABC_ARabc[1:Sp,1]),
  Distribution = factor(c(rep("EasyABC", Sp), rep("ABC", Sp), rep("ABCRegAd", Sp), rep("ABCabc", Sp))))

dentheta <- ggplot(df1, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) + 
  geom_vline(xintercept = theta1, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = TeX("Posterior density plot: $theta$"), x = TeX("$theta$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red", "green", "purple")) +  theme_minimal() +
  theme(legend.title = element_blank())

df2 <- data.frame(
  Value = c(PostABC_AR[1:Sp,4], PostABC_ARown[1:Sp,4], PostABC_ARownRegAd[1:Sp,4], PostABC_ARabc[1:Sp,4]),
  Distribution = factor(c(rep("EasyABC", Sp), rep("Own", Sp), rep("ABCRegAd", Sp), rep("ABCabc", Sp))))

deng <- ggplot(df2, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) + 
  geom_vline(xintercept = g, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = TeX("Posterior density plot: g"), x = TeX("$g$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red", "green", "purple")) +  theme_minimal() +
  theme(legend.title = element_blank())

df3 <- data.frame(
  Value = c(PostABC_AR[1:Sp,5], PostABC_ARown[1:Sp,5], PostABC_ARownRegAd[1:Sp,5], PostABC_ARabc[1:Sp,4]),
  Distribution = factor(c(rep("EasyABC", Sp), rep("Own", Sp), rep("ABCRegAd", Sp), rep("ABCabc", Sp))))

denk <- ggplot(df3, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +
  geom_vline(xintercept = k, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = TeX("Posterior density plot: k"), x = TeX("$k$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red", "green", "purple")) +  theme_minimal() +
  theme(legend.title = element_blank())

library(ggpubr)
ggarrange(dentheta, deng, denk, labels = c("A", "B", "C"), ncol = 3, nrow = 1,
          legend = "bottom", common.legend = TRUE)

######## BSL: Financial returns application ############
rm(list = ls()); set.seed(010101)
library(BSL)
dfExcRate <- read.csv(file = "https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/ExchangeRate.csv", sep = ",", header = T)
attach(dfExcRate)
str(dfExcRate)
y <- USDEUR
n <- length(y)
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
par0 <- c(0.2, 1, 1, 0, 0.8) # c(0.5, 2, 1, 0, 1) 
Modelgk <- newModel(fnSim = RGKnew, fnSum = SumSt, theta0 = par0, fnLogPrior = LogPrior, verbose = FALSE)
validObject(Modelgk)
M <- 200 # 500 Number of iterations to calculate mu and sigma
S <- 1100 # 11000 Number of MCMC iterations
burnin <- 100 # 1000 # Burn in iterations
thin <- 2 # 10 # Thining parameter
keep <- seq(burnin + 1, S, thin)
tune <- 0.01 # Tuning parameter RW MH
simgk <- simulation(Modelgk, n = M, theta = par0, seed = 10)
par(mfrow = c(4, 3))
# Check if the summary statistics are roughly normal
for (i in 1:12){
  eval <- seq(min(simgk$ssx[, i]), max(simgk$ssx[, i]), 0.001)
  densnorm <- dnorm(eval, mean = mean(simgk$ssx[, i]), sd(simgk$ssx[, i])) 
  plot(density(simgk$ssx[, i]), main = "", xlab = "")
  lines(eval, densnorm, col = "red")
} 
Lims <- matrix(c(-1, 0, 0, -5, -0.5, 1, rep(5, 4)), 5, 2)
tick <- Sys.time()
Resultsgk <- bsl(y = y, n = M, M = S, model = Modelgk, covRandWalk = tune*diag(5),
                 method = "BSL", thetaNames = expression(theta, a, b, g, k), 
                 logitTransformBound = Lims, plotOnTheFly = TRUE)
tock <- Sys.time()
tock - tick
PostChain <- coda::mcmc(Resultsgk@theta[keep,])
summary(PostChain)
CovarRWnew <- var(PostChain)
M <- 500 # Number of iterations to calculate mu and sigma
S <- 30000 # 11000 # 11000 Number of MCMC iterations
burnin <- 10000 # 1000 # Burn in iterations
thin <- 10 # 10 # Thining parameter
keep <- seq(burnin + 1, S, thin)
tune <- 1 # Tuning parameter RW MH
get_mode <- function(x) {
  uniq_vals <- unique(x)
  uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
}
# apply(PostChain, 2, get_mode)
ModelgkNew <- newModel(fnSim = RGKnew, fnSum = SumSt, theta0 = par0, fnLogPrior = LogPrior, verbose = FALSE)
logitTransform <- function(par, a, b){
  logtrans <- log((par - a)/(b - par))
  return(logtrans)
}
ParTrans <- matrix(NA, dim(PostChain)[1], 5)
for(j in 1:5){
  ParTrans[,j] <- logitTransform(par = PostChain[,j], a = Lims[j,1], b = Lims[j,2])
}
CovarRW <- var(ParTrans)
# The acceptance rate is really low in this setting!!!
# tick <- Sys.time()
# ResultsgkNew <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
#                     method = "BSL", thetaNames = expression(theta, a, b, g, k), 
#                     logitTransformBound = Lims, plotOnTheFly = TRUE)
# tock <- Sys.time()
# tock - tick
# PostChain <- coda::mcmc(ResultsgkNew@theta[keep,])
# plot(PostChain)
# ResultsgkNew@acceptanceRate
# plot(ResultsgkNew@loglike[keep], type = "l")
# sd(ResultsgkNew@loglike[keep])

# These two didn't work in this example!!!
# tick <- Sys.time()
# ResultsgkuBLS <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
#                     method = "uBSL", thetaNames = expression(theta, a, b, g, k), 
#                     logitTransformBound = Lims, plotOnTheFly = TRUE)
# tock <- Sys.time()
# tock - tick
# PostChainuBLS <- coda::mcmc(ResultsgkuBLS@theta[keep,])
# plot(PostChainuBLS)
# ResultsgkuBLS@acceptanceRate
# plot(ResultsgkuBLS@loglike[keep], type = "l")
# sd(ResultsgkuBLS@loglike[keep])
# 
# tick <- Sys.time()
# ResultsgksemiBLS <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
#                      method = "semiBSL", thetaNames = expression(theta, a, b, g, k), 
#                      logitTransformBound = Lims, plotOnTheFly = TRUE)
# tock <- Sys.time()
# tock - tick
# PostChainsemiBLS <- coda::mcmc(ResultsgksemiBLS@theta[keep,])
# plot(PostChainsemiBLS)
# ResultsgksemiBLS@acceptanceRate
# plot(ResultsgksemiBLS@loglike[keep], type = "l")
# sd(ResultsgksemiBLS@loglike[keep])

# The acceptance rate in this is really low!!!
# tick <- Sys.time()
# ResultsgkmisspecBLS <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
#                         method = "BSLmisspec", thetaNames = expression(theta, a, b, g, k), 
#                         logitTransformBound = Lims, misspecType = "mean", tau = 0.5, plotOnTheFly = TRUE)
# tock <- Sys.time()
# tock - tick
# PostChainmisspecBLS <- coda::mcmc(ResultsgkmisspecBLS@theta[keep,])
# plot(PostChainmisspecBLS)
# ResultsgkmisspecBLS@acceptanceRate
# plot(ResultsgkmisspecBLS@loglike[keep], type = "l")
# sd(ResultsgkmisspecBLS@loglike[keep])

tick <- Sys.time()
ResultsgkmisspecBLSvar <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
                           method = "BSLmisspec", thetaNames = expression(theta, a, b, g, k), 
                           logitTransformBound = Lims, misspecType = "variance", tau = 0.5, plotOnTheFly = TRUE)
tock <- Sys.time()
tock - tick
PostChainmisspecBLSvar <- coda::mcmc(ResultsgkmisspecBLSvar@theta[keep,])
# plot(PostChainmisspecBLSvar)
# ResultsgkmisspecBLSvar@acceptanceRate
# plot(ResultsgkmisspecBLSvar@loglike[keep], type = "l")
# sd(ResultsgkmisspecBLSvar@loglike[keep])
# Figures 
library(ggplot2); library(latex2exp)
Sp <- length(keep)
df1 <- data.frame(
  Value = c(PostChainmisspecBLSvar[,1]),
  Distribution = factor(c(rep("BSL", Sp)))
)

dentheta <- ggplot(df1, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = TeX("Posterior density plot: $theta$"), x = TeX("$theta$"), y = "Posterior density") +
  scale_color_manual(values = c("blue")) +  theme_minimal() +
  theme(legend.title = element_blank())

df2 <- data.frame(
  Value = c(PostChainmisspecBLSvar[,4]),
  Distribution = factor(c(rep("BSL", Sp)))
)

deng <- ggplot(df2, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = "Posterior density plot: g", x = "g", y = "Posterior density") +
  scale_color_manual(values = c("blue")) +  theme_minimal() +
  theme(legend.title = element_blank())

df3 <- data.frame(
  Value = c(PostChainmisspecBLSvar[,5]),
  Distribution = factor(c(rep("BSL", Sp)))
)

denk <- ggplot(df3, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = "Posterior density plot: k", x = "k", y = "Posterior density") +
  scale_color_manual(values = c("blue")) +  theme_minimal() +
  theme(legend.title = element_blank())

library(ggpubr)
ggarrange(dentheta, deng, denk, labels = c("A", "B", "C"), ncol = 3, nrow = 1,
          legend = "bottom", common.legend = TRUE)


#####################################################################
tick <- Sys.time()
ResultsgkmisspecBLSvarnew <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRWnew,
                              method = "BSLmisspec", thetaNames = expression(theta, a, b, g, k), 
                              misspecType = "variance", tau = 0.5, plotOnTheFly = TRUE)
tock <- Sys.time()
tock - tick
keep <- seq(burnin + 1, S, thin) # keep <- seq(burnin + 1, S, thin)
PostChainmisspecBLSvarnew <- coda::mcmc(ResultsgkmisspecBLSvarnew@theta[keep,])
# plot(PostChainmisspecBLSvarnew)
# ResultsgkmisspecBLSvarnew@acceptanceRate
# plot(ResultsgkmisspecBLSvarnew@loglike[keep], type = "l")
# sd(ResultsgkmisspecBLSvarnew@loglike[keep])

# Figures 
library(ggplot2); library(latex2exp)
Sp <- length(keep)
df1 <- data.frame(
  Value = c(PostChainmisspecBLSvarnew[,1]),
  Distribution = factor(c(rep("BSL", Sp)))
)

dentheta <- ggplot(df1, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = TeX("Posterior density plot: $theta$"), x = TeX("$theta$"), y = "Posterior density") +
  scale_color_manual(values = c("blue")) +  theme_minimal() +
  theme(legend.title = element_blank())

df2 <- data.frame(
  Value = c(PostChainmisspecBLSvarnew[,4]),
  Distribution = factor(c(rep("BSL", Sp)))
)

deng <- ggplot(df2, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = "Posterior density plot: g", x = "g", y = "Posterior density") +
  scale_color_manual(values = c("blue")) +  theme_minimal() +
  theme(legend.title = element_blank())

df3 <- data.frame(
  Value = c(PostChainmisspecBLSvarnew[,5]),
  Distribution = factor(c(rep("BSL", Sp)))
)

denk <- ggplot(df3, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = "Posterior density plot: k", x = "k", y = "Posterior density") +
  scale_color_manual(values = c("blue")) +  theme_minimal() +
  theme(legend.title = element_blank())

library(ggpubr)
ggarrange(dentheta, deng, denk, labels = c("A", "B", "C"), ncol = 3, nrow = 1,
          legend = "bottom", common.legend = TRUE)
