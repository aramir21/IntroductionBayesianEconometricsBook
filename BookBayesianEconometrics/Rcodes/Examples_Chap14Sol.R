# g-and-k distribution
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
theta1 <- 0.6; a <- 1; b <- 0.5; g <- -1; k <- 1
parpop <- c(theta1, a, b, g, k)
n <- 250
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
tick <- Sys.time()
ResABC <- ABC(S = S, a = 0.001, y = y)
tock <- Sys.time()
tock - tick
PostABC_ARown <- ResABC[["SelPrior"]]
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
                        summary_stat_target = sum_stat_obs, nb_simul=100000, tol = 0.01,
                        progress_bar = TRUE)
tock <- Sys.time()
tock - tick
PostABC_AR <- coda::mcmc(ABC_AR$param)

# Figures 
library(ggplot2); library(latex2exp)
df1 <- data.frame(
  Value = c(PostABC_AR[,1], PostABC_ARown[,1]),
  Distribution = factor(c(rep("EasyABC", 1000), rep("Own", 1000))))

dentheta <- ggplot(df1, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = TeX("Density plot posterior: $theta$"), x = TeX("$theta$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red")) +  theme_minimal() +
  theme(legend.title = element_blank())

df2 <- data.frame(
  Value = c(PostABC_AR[,4], PostABC_ARown[,4]),
  Distribution = factor(c(rep("EasyABC", 1000), rep("Own", 1000))))

deng <- ggplot(df2, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = TeX("Density plot posterior: g"), x = TeX("$g$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red")) +  theme_minimal() +
  theme(legend.title = element_blank())

df3 <- data.frame(
  Value = c(PostABC_AR[,5], PostABC_ARown[,5]),
  Distribution = factor(c(rep("EasyABC", 1000), rep("Own", 1000))))

denk <- ggplot(df3, aes(x = Value, color = Distribution)) +   geom_density(linewidth = 1) +  
  labs(title = TeX("Density plot posterior: k"), x = TeX("$k$"), y = "Posterior density") +
  scale_color_manual(values = c("blue", "red")) +  theme_minimal() +
  theme(legend.title = element_blank())

library(ggpubr)
ggarrange(dentheta, deng, denk, labels = c("A", "B", "C"), ncol = 3, nrow = 1,
          legend = "bottom", common.legend = TRUE)

