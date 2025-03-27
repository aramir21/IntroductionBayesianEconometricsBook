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
