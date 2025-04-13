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
  Etay <- c(Oct[4], eta1, eta2, eta3, autocor[["acf"]][2:3])
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
S <- 1100 # Number of MCMC iterations
burnin <- 100 # Burn in iterations
thin <- 2 # Thining parameter
keep <- seq(burnin + 1, S, thin)
par0 <- c(0.5, 2, 1, 0, 1) 
Modelgk <- newModel(fnSim = RGKnew, fnSum = SumSt, theta0 = par0, fnLogPrior = LogPrior, verbose = FALSE)
validObject(Modelgk)
keep <- seq(burnin + 1, S, thin)
tune <- 0.1 # Tuning parameter RW MH
simgk <- simulation(Modelgk, n = M, theta = par0, seed = 10)
par(mfrow = c(2, 3))
# Check if the summary statistics are roughly normal
for (i in 1:6){
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
CovarRWnew <- var(PostChain)
M <- 500 # Number of iterations to calculate mu and sigma
S <- 6000 # Number of MCMC iterations
burnin <- 1000 # Burn in iterations
thin <- 5 # Thining parameter
keep <- seq(burnin + 1, S, thin)
tune <- 1 # Tuning parameter RW MH
# get_mode <- function(x) {
#   uniq_vals <- unique(x)
#   uniq_vals[which.max(tabulate(match(x, uniq_vals)))]
# }
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
# Vanilla BLS
tick <- Sys.time()
ResultsgkVanila <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
                       method = "BSL", thetaNames = expression(theta, a, b, g, k),
                       logitTransformBound = Lims, plotOnTheFly = TRUE)
tock <- Sys.time()
tock - tick
# Unbiased BLS
tick <- Sys.time()
ResultsgkuBLS <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
                     method = "uBSL", thetaNames = expression(theta, a, b, g, k),
                     logitTransformBound = Lims, plotOnTheFly = TRUE)
tock <- Sys.time()
tock - tick
# semi-parametric BLS
tick <- Sys.time()
ResultsgksemiBLS <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
                        method = "semiBSL", thetaNames = expression(theta, a, b, g, k),
                        logitTransformBound = Lims, plotOnTheFly = TRUE)
tock <- Sys.time()
tock - tick
# Misspecified BLS in mean
tick <- Sys.time()
ResultsgkmisspecBLSmean <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
                               method = "BSLmisspec", thetaNames = expression(theta, a, b, g, k),
                               logitTransformBound = Lims, misspecType = "mean", tau = 0.5, plotOnTheFly = TRUE)
tock <- Sys.time()
tock - tick
# Misspecified BLS in variance
tick <- Sys.time()
ResultsgkmisspecBLSvar <- bsl(y = y, n = M, M = S, model = ModelgkNew, covRandWalk = tune*CovarRW,
                              method = "BSLmisspec", thetaNames = expression(theta, a, b, g, k), 
                              logitTransformBound = Lims, misspecType = "variance", tau = 0.5, plotOnTheFly = TRUE)
tock <- Sys.time()
tock - tick
# Summary all models
gkResults <- list(ResultsgkVanila, ResultsgkuBLS, ResultsgkmisspecBLSmean, ResultsgkmisspecBLSvar)
names(gkResults) <- c("BSL", "uBSL", "missBSLmean", "missBSLvar")
t(sapply(gkResults, summary))
combinePlotsBSL(gkResults, which = 1, thin = thin, burnin = burnin)
save(gkResults, file = "ExerBSLgk.RData")
# load(file = "ExerBSLgk.RData")
combinePlotsBSL(gkResults, which = 2, thetaTrue = parpop,
                thin = thin, options.linetype = list(values = 1:8),
                options.size = list(values = rep(1, 8)),
                options.theme = list(plot.margin = grid::unit(rep(0.03, 4), "npc"),
                axis.title = ggplot2::element_text(size = 12),
                axis.text = ggplot2::element_text(size = 8),
                legend.text = ggplot2::element_text(size = 12)))
# semi BSL and miss BSL mean are totally diferent to other methods 
gkResultsNew <- list(gkResults[["BSL"]], gkResults[["uBSL"]], gkResults[["missBSLvar"]])
names(gkResultsNew) <- c("BSL", "uBSL", "missBSLvar")
combinePlotsBSL(gkResultsNew, which = 2, thetaTrue = parpop,
                thin = thin, options.linetype = list(values = 1:8),
                options.size = list(values = rep(1, 8)),
                options.theme = list(plot.margin = grid::unit(rep(0.03, 4), "npc"),
                legend.text = ggplot2::element_text(size = 12)))
# All 95% credible intervals encompass the population parameters.

########################## Simulation exercise: Linear regression SVI using Stan ########################
library(rstan)
# --- Simulated data ---
N <- 100000
K <- 2
B <- rep(1, K + 1)
X <- cbind(1, matrix(rnorm(N * K), N, K))
sig2 <- 1
y <- X %*% B + rnorm(N, 0, sig2^0.5)
# --- Priors ---
b0 <- 0; B0 <- 1000; a0 <- 0.01; d0 <- 0.01
# --- Prepare data for Stan ---
data_list <- list(
  N = N,
  K = K + 1,  # K+1 because of the intercept
  X = X,
  y = as.vector(y),
  b0 = b0,
  B0 = B0,
  a0 = a0,
  d0 = d0
)
# --- Stan model compilation ---
stan_model_code <- "
data {
  int<lower=0> N;             // number of data points
  int<lower=0> K;             // number of predictors
  matrix[N, K] X;             // predictor matrix
  vector[N] y;                // response vector
  real b0;                    // location parameter for the Normal prior on beta
  real<lower=0> B0;           // variance parameter for the Normal prior on beta
  real<lower=0> a0;           // shape parameter for the Inv-Gamma prior on sigma^2
  real<lower=0> d0;           // rate parameter for the Inv-Gamma prior on sigma^2
}

parameters {
  vector[K] beta;             // coefficients of the linear model
  real<lower=0> sigma2;       // variance parameter (sigma^2)
}

model {
  // Priors
  beta ~ normal(b0, B0);       // Normal prior on beta with variance 100
  sigma2 ~ inv_gamma(a0/2, d0/2); // Inverse Gamma prior on sigma^2

  // Likelihood
  y ~ normal(X * beta, sqrt(sigma2)); // Normal likelihood with linear predictor
}

generated quantities {
  vector[N] y_pred;

  for (n in 1:N) {
    y_pred[n] = normal_rng(dot_product(X[n], beta), sqrt(sigma2));
  }
}
"
# Compile Stan model
stan_model <- stan_model(model_code = stan_model_code)
# --- Perform Stochastic Variational Inference (SVI) ---
fit <- vb(stan_model, data = data_list, algorithm = "meanfield", iter = 10000)
# --- Extract the results ---
mu_estimate <- extract(fit)$beta
summary(coda::mcmc(mu_estimate))
sigma2_estimate <- extract(fit)$sigma2
summary(coda::mcmc(as.vector(sigma2_estimate)))

########################## Simulation exercise: Gaussian mixture: 2 components ########################## 
rm(list = ls())
set.seed(010101)
# Simulate data from a 2-component mixture model
n <- 1000
K <- 2; H <- 2
x <- rnorm(n)
z <- rbinom(n, 1, 0.5)  # Latent class indicator
y <- ifelse(z == 0, rnorm(n, 2 + 1.5*x, 1), rnorm(n, -1 + 0.5*x, 0.8))
X <- cbind(1,x)
# Iterations
S <- 100
##### From scratch ######
# Hyperparameters
b0h <- rep(0, K); B0h <- diag(1000,K); B0hi <- solve(B0h) 
a0h <- 0.01; d0h <- 0.01; a0 <- rep(1/H, H)
# Innitial values
dhn <- c(50, 50); ahn <- c(500, 500)
bhn <- matrix(c(2, 1.5, -1, 0.5), K, H)
Bhn <- array(c(1,0,0,1,1,0,0,1), c(K, K, H))
an <- c(500, 500)
VarParPsi <- function(dhn, ahn, bhn, Bhn, an){
  logrhoih <- matrix(NA, n, H)
  for(i in 1:n){
    for(h in 1:H){
      logrhoih[i,h] <- -0.5*log(2*pi) - 0.5*(log(dhn[h]/2) - digamma(ahn[h]/2)) -
        0.5*(ahn[h]/dhn[h])*((y[i] - t(X[i,])%*%bhn[,h])^2 + t(X[i,])%*%Bhn[,,h]%*%X[i,]) +
        digamma(an[h]) - digamma(sum(an))
    }
  }
  rhoih <- exp(logrhoih - apply(logrhoih, 1, max)) 
  rih <- rhoih / rowSums(rhoih)
  return(rih)
}
Rih <- VarParPsi(dhn = dhn, ahn = ahn, bhn = bhn, Bhn = Bhn, an = an)
VarParLambda <- function(rih){
  an <- a0 + colSums(rih)
  return(an)
}
VarParLambda(rih = Rih)
VarParBetah <- function(ahn, dhn, rh){
  rXX <- matrix(0, K, K)
  rXy <- matrix(0, K, 1)
  for(i in 1:n){
    rXXi <- rh[i]*X[i,]%*%t(X[i,])
    rXX <- rXX + rXXi
    rXyi <- rh[i]*X[i,]*y[i]
    rXy <- rXy + rXyi
  }
  # Bhn <- solve(B0hi + ahn/dhn * rXX + diag(1e-8, K))
  Bhn <- solve(B0hi + ahn/dhn * rXX)
  bhn <- Bhn%*%(B0hi%*%b0h + ahn/dhn * rXy)
  return(list(bhn = bhn, Bhn = Bhn))
}
h <- 2
VarParBetah(ahn = ahn[h], dhn = dhn[h], rh = Rih[,h])
VarParSigma2h <- function(bhn, Bhn, rh){
  ahn <- a0h + sum(rh)
  dhnterm <- 0
  for(i in 1:n){
    dhntermi <- rh[i]*((y[i] - t(X[i,])%*%bhn)^2 + t(X[i,])%*%Bhn%*%X[i,])
    dhnterm <- dhnterm + dhntermi
  }
  dhn <- dhnterm + d0h
  return(list(ahn = ahn, dhn = dhn))
}
h <- 1
VarParSigma2h(bhn = bhn[,h], Bhn = Bhn[,,h], rh = Rih[,h])
# Auxiliar function to calculate the log multivariate beta function
# Define Gamma and Beta functions
multivariate_beta <- function(alpha) {
  prod_gammas <- sum(sapply(alpha, lgamma))
  sum_alpha <- sum(alpha)
  gamma_sum_alpha <- lgamma(sum_alpha)
  return(prod_gammas - gamma_sum_alpha)
}
ELBO <- function(rih, dhn, ahn, bhn, Bhn, an){
  LogLikh <- 0; LogPrioPsih <- 0
  LogPriorBeta <- 0; LogPriorSig2 <- 0
  for(h in 1:H){
    for(i in 1:n){
      LogLikih <- rih[i, h]*(-0.5*log(2*pi) - 0.5*(log(dhn[h]/2) - digamma(ahn[h]/2)) -
                               0.5*(ahn[h]/dhn[h])*((y[i] - t(X[i,])%*%bhn[,h])^2 + 
                                                      t(X[i,])%*%Bhn[,,h]%*%X[i,]))
      LogLikh <- LogLikh + LogLikih
      LogPrioPsiih <- rih[i, h]*(digamma(an[h]) - digamma(sum(an)))
      LogPrioPsih <- LogPrioPsih + LogPrioPsiih
    }
    LogPriorBetah <- -0.5*K*log(2*pi) - 0.5*log(det(B0h)) - sum(diag(Bhn[,,h]%*%B0hi)) -
      t(bhn[,h] - b0h)%*%B0hi%*%(bhn[,h] - b0h)
    LogPriorBeta <- LogPriorBeta + LogPriorBetah
    LogPriorSig2h <- -(a0h/2+1)*(log(dhn[h]/2) - digamma(ahn[h]/2)) - 0.5*d0h*(ahn[h]/dhn[h]) +
      (a0h/2)*log(d0h/2) - lgamma(a0h/2)
    LogPriorSig2 <- LogPriorSig2 + LogPriorSig2h  
  }
  LogJoint <- LogLikh + LogPrioPsih + LogPriorBeta + LogPriorSig2
  logqPsi <- 0; logqLambda <- 0 
  logqBeta <- 0; logqSig2 <- 0
  for(h in 1:H){
    for(i in 1:n){
      logqPsiih <- rih[i, h] * log(rih[i, h])
      logqPsi <- logqPsi + logqPsiih
    }
    logqLambdah <- (an[h] - 1)*(digamma(an[h]) - digamma(sum(an)))
    logqLambda <- logqLambda + logqLambdah
    logqBetah <- -0.5*K*log(2*pi + 1) - 0.5*log(det(Bhn[,,h])) 
    logqBeta <- logqBeta + logqBetah
    logqSig2h <- -log(dhn[h]/2) + (ahn[h]/2+1)*digamma(ahn[h]/2) - ahn[h]/2 -
      lgamma(ahn[h]/2)
    logqSig2 <- logqSig2 + logqSig2h
  }
  logqLambda <- logqLambda - multivariate_beta(alpha = an)
  logq <- logqPsi + logqLambda + logqBeta + logqSig2
  ELBO <- LogJoint - logq
  return(ELBO)
}
ELBO(rih = Rih, dhn = dhn, ahn = ahn, bhn = bhn, Bhn = Bhn, an = an)
# Initial values
# Use rough values from data (not cheating, just guiding)
fit_lm1 <- lm(y ~ x)
bhn[,1] <- coef(fit_lm1)
bhn[,2] <- coef(fit_lm1) + rnorm(K, 0, 0.5)  # small perturbation

Bhn[,,1] <- diag(0.1, K)
Bhn[,,2] <- diag(0.1, K)

ahn <- rep(10, H)
dhn <- rep(10, H)

an <- c(300, 700)  # break symmetry (prior pseudo-counts)
ELBOs <- rep(NA, S)
ELBOs[1] <- ELBO(rih = Rih, dhn = dhn, ahn = ahn, bhn = bhn, Bhn = Bhn, an = an)
epsilon <- 1e-5
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
tick <- Sys.time()
for(s in 2:S){
  Rih <- VarParPsi(dhn = dhn, ahn = ahn, bhn = bhn, Bhn = Bhn, an = an)
  an <- VarParLambda(rih = Rih)
  bhn <- matrix(NA, K, H)
  Bhn <- array(NA, c(K, K, H))
  for(h in 1:H){
    EnsVarParBetah <- VarParBetah(ahn = ahn[h], dhn = dhn[h], rh = Rih[,h])
    bhn[,h] <- EnsVarParBetah[["bhn"]]
    Bhn[,,h] <- EnsVarParBetah[["Bhn"]] 
  }
  ahn <- rep(NA, H)
  dhn <- rep(NA, H)
  for(h in 1:H){
    EnsVarParSig2h <- VarParSigma2h(bhn = bhn[,h], Bhn = Bhn[,,h], rh = Rih[,h])
    ahn[h] <- EnsVarParSig2h[["ahn"]]
    dhn[h] <- EnsVarParSig2h[["dhn"]] 
  }
  ELBOs[s] <- ELBO(rih = Rih, dhn = dhn, ahn = ahn, bhn = bhn, Bhn = Bhn, an = an)
  # Check if lower bound decreases
  if (ELBOs[s] < ELBOs[s - 1]) { message("Lower bound decreases!\n")}
  # Check for convergence
  if (ELBOs[s] - ELBOs[s - 1] < epsilon) { break }
  # Check if VB converged in the given maximum iterations
  if (s == S) {warning("VB did not converge!\n")}
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0), "% done"))
}
close(pb)
tock <- Sys.time()
tock - tick
plot(ELBOs[1:s], type = "l", main = "ELBO over iterations", ylab = "ELBO", xlab = "Iteration")
Sig21 <- MCMCpack::rinvgamma(S, ahn[1]/2, dhn[1]/2) 
summary(coda::mcmc(Sig21))
Sig22 <- MCMCpack::rinvgamma(S, ahn[2]/2, dhn[2]/2) 
summary(coda::mcmc(Sig22))
Beta1 <- MASS::mvrnorm(S, mu = bhn[,1], Sigma = Bhn[,,1])
summary(coda::mcmc(Beta1))
Beta2 <- MASS::mvrnorm(S, mu = bhn[,2], Sigma = Bhn[,,2])
summary(coda::mcmc(Beta2))
Lambda <- MCMCpack::rdirichlet(S, alpha = an)
summary(coda::mcmc(Lambda))
