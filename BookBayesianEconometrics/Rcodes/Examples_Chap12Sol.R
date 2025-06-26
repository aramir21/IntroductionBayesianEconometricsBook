######### Bayesian LASSO with hierarchical lambda #########
rm(list = ls()); set.seed(10101)
library(monomvn)

# Parameters
n <- 500  # sample size
k <- 100  # number of predictors
s <- 10   # number of non-zero coefficients

# Generate design matrix
X <- matrix(rnorm(n * k), nrow = n, ncol = k)

# True beta: first s coefficients are non-zero, rest are zero
beta_true <- c(runif(s, -3, 3), rep(0, k - s))

# Generate response with some noise
sigma <- 1
y <- X %*% beta_true + rnorm(n, sd = sigma)
df <- data.frame(X,y)

### Using monomvn ###
library(monomvn)
library(statmod)

# Fit model using Bayesian LASSO
fit <- monomvn::blasso(X, y, T = 5000, verb = 0)

# burn-in removal
burnin <- 1000
beta_samples <- fit$beta[(burnin + 1):5000, ]  # exclude burn-in
intercept_samples <- fit$mu[(burnin + 1):5000] # intercept posterior draws

# Posterior mean of coefficients
beta_post_mean <- colMeans(beta_samples)

# Intercept summary
mean(intercept_samples)
quantile(intercept_samples, c(0.025, 0.975))

# Local shrinkage parameters
colMeans(is.na(fit$tau2i))  # proportion of NA values per coefficient

#### Gibss sampler ####
# Hyperparameters
r0 <- 1
d0 <- 1
# MCMC parameters
mcmc <- 5000
burnin <- 5000
tot <- mcmc + burnin
thin <- 1
# Posterior distributions programming the Gibbs sampling
# Standardized variables
W <- scale(X)
yh <- y - mean(y)
# Gibbs sampling functions
PostBeta <- function(sig2, tau2){
  Bn <- solve(diag(1/tau2) + t(W)%*%W)
  bn <- Bn%*%t(W)%*%yh
  Beta <- MASS::mvrnorm(1, bn, sig2*Bn)
  return(Beta)
}

PostSig2 <- function(Beta, tau2){
  an <- n - 1 + k
  dn <- t(yh - W%*%Beta)%*%(yh - W%*%Beta) + t(Beta)%*%diag(1/tau2)%*%Beta
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}

PostTau2 <- function(sig2, Beta, lambda, p){
  mukn <- (lambda^2*sig2/Beta[p]^2)^0.5
  lambdan <- lambda^2
  tauk2k <- rinvgauss(1, mean = mukn, shape = lambdan)
  return(tauk2k)
}

PostLambda2 <- function(tau2){
  sh <- k + r0
  rat <- 0.5 * sum(tau2) + d0
  lambda2 <- rgamma(1, shape = sh, rate = rat)
  return(lambda2^0.5)
}

PostBeta0 <- function(sig2){
  beta0 <- rnorm(1, mean = mean(y), sd = (sig2/n)^0.5)
  return(beta0)
}

PostBetas <- matrix(0, mcmc+burnin, k)
PostSigma2 <- rep(0, mcmc+burnin)
PostTaus2 <- matrix(0, mcmc+burnin, k)
PostLambda <- rep(0, mcmc+burnin)
PostBetas0 <- rep(0, mcmc+burnin)
Beta <- rep(0, k)
tau2 <- rep(1, k)
lambda <- 1
# create progress bar in case that you want to see iterations progress
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  sig2 <- PostSig2(Beta = Beta, tau2 = tau2)
  PostSigma2[s] <- sig2
  Beta <- PostBeta(sig2 = sig2, tau2 = tau2)
  PostBetas[s,] <- Beta
  tau2 <- sapply(1:k, function(i){PostTau2(sig2 = sig2, Beta = Beta, lambda = lambda, p = i)})
  PostTaus2[s,] <- tau2
  lambda <- PostLambda2(tau2 = tau2)
  PostLambda[s] <- lambda
  beta0 <- PostBeta0(sig2 = sig2)
  PostBetas0[s] <- beta0
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0), "% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
summary(coda::mcmc(PosteriorBetas))
PosteriorBeta0 <- PostBetas0[keep]
summary(coda::mcmc(PosteriorBeta0))
PosteriorSigma2 <- PostSigma2[keep]
summary(coda::mcmc(PosteriorSigma2))
PosteriorTau2s <- PostTaus2[keep,]
summary(coda::mcmc(PosteriorTau2s))
PosteriorLambda <- PostLambda[keep]
summary(coda::mcmc(PosteriorLambda))

# Posterior mean of coefficients
beta_post_meanNew <- colMeans(PosteriorBetas)

# Plot true vs estimated betas
plot(beta_true, beta_post_mean, pch = 19, col = "steelblue",
     xlab = "True beta", ylab = "Posterior mean of beta",
     main = "Bayesian LASSO (blasso): Shrinkage Effect")
abline(0, 1, col = "red", lty = 2)

# Add another set of posterior means
points(beta_true, beta_post_meanNew, pch = 17, col = "darkorange")

# Add legend
legend("topleft", legend = c("blasso", "Our Gibbs"), 
       col = c("steelblue", "darkorange"), 
       pch = c(19, 17), bty = "n")


######### Application SSVS: Conflict #########
rm(list = ls())
set.seed(010101)
Data <- read.csv("https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/Conflict.csv", sep = ",", header = TRUE, quote = "")

# Scale regressors
W <- as.matrix(scale(Data[, -c(1, 2)]))
niter <- 15000
y <- unlist(Data[, 2])

# Linear model
SSBoomLinear <- lm.spike(y ~ W, niter = niter)
Models <- SSBoomLinear$beta != 0
PIP <- colMeans(SSBoomLinear$beta != 0)
# Convert the logical matrix to a data frame and then to a tibble
df <- as.data.frame(Models); df_tbl <- as_tibble(df)
# Count identical rows
row_counts <- df_tbl %>% count(across(everything()), name = "frequency") %>% arrange(desc(frequency))
colnames(W)
# Get indices of best model
for(l in 1:5) {
  print(which(row_counts[l,] == 1))
  print(row_counts[l, dim(row_counts)[2]]/niter)
  }
# Coefficients
SummarySS <- summary(coda::mcmc(SSBoomLinear$beta))
SummarySS

# Logit model
SSBoomLogit <- logit.spike(y ~ W, niter = niter)
ModelsLogit <- SSBoomLogit$beta != 0
PIPLogit <- colMeans(SSBoomLogit$beta != 0)
# Convert the logical matrix to a data frame and then to a tibble
df <- as.data.frame(ModelsLogit); df_tbl <- as_tibble(df)
# Count identical rows
row_counts <- df_tbl %>% count(across(everything()), name = "frequency") %>% arrange(desc(frequency))
colnames(W)
# Get indices of best model
for(l in 1:5) {
  print(which(row_counts[l,] == 1))
  print(row_counts[l, dim(row_counts)[2]]/niter)
}
# Coefficients
SummarySSlogit <- summary(coda::mcmc(SSBoomLogit$beta))
SummarySSlogit

######### BLASSO and SSVS: K > N #########
rm(list = ls()); set.seed(10101)
library(BoomSpikeSlab)
library(dplyr)
library(tibble)
# Parameters
n <- 500  # sample size
k <- 600  # number of predictors
s <- 10   # number of non-zero coefficients
# Generate design matrix
X <- matrix(rnorm(n * k), nrow = n, ncol = k)
# True beta: first s coefficients are non-zero, rest are zero
beta_true <- c(runif(s, -3, 3), rep(0, k - s))
# Generate response with some noise
sigma <- 1
y <- X %*% beta_true + rnorm(n, sd = sigma)
df <- data.frame(X,y)
### Using BoomSpikeSlab ###
#Scale regressors
W <- scale(X); yh <- y - mean(y)
niter <- 5000
######Estimate model########
# Least squared
Reg <- lm(y ~ X)
summary(Reg)
## BLASSO ##
# Fit the model
fit <- bayesreg::bayesreg(y ~ X, data = df, model = "gaussian", prior = "lasso", 
                          n.samples = niter, burnin = 1000)
# Check summary
summary(fit)
# Extract posterior means of beta
beta_post_meanLASSO <- rowMeans(fit$beta)

## SSVS ##
SSBoomNew <- lm.spike(yh ~ W - 1, niter = niter)
# Coefficients
SummarySS <- summary(coda::mcmc(SSBoomNew$beta))
# Extract posterior means of beta
beta_post_meanSSVS <- SummarySS$statistics[, 1]

# Compare true vs estimated
plot(beta_true, beta_post_meanLASSO, pch = 19, col = "steelblue",
     xlab = "True beta", ylab = "Posterior mean of beta",
     main = "Bayesian LASSO and SSVS: K > N")
abline(0, 1, col = "red", lty = 2)
# Add another set of posterior means
points(beta_true, beta_post_meanSSVS, pch = 17, col = "darkorange")
# Add legend
legend("topleft", legend = c("BLASSO", "SSVS"), 
       col = c("steelblue", "darkorange"), 
       pch = c(19, 17), bty = "n")

cbind(beta_true, beta_post_meanLASSO, beta_post_meanSSVS)

####### BART #########
rm(list = ls()); set.seed(10101)
library(BART); library(ggplot2); library(dplyr); library(tidyr)
N <- 500; K <- 10
# Simulate the dataset
MeanFunct <- function(x){
  f <- 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
  return(f)
}
sig2 <- 1
e <- rnorm(N, 0, sig2^0.5)
X <- matrix(runif(N*K),N,K)
y <- MeanFunct(X[,1:5]) + e

# Hyperparameters
alpha <- 0.95; beta <- 2; k <- 2
v <- 3; q <- 0.9; J <- 200
# MCMC parameters
MCMCiter <- 1000; burnin <- 100; thinning <- 1


# Partial dependence function
BARTfit <- wbart(x.train = X, y.train = y, base = alpha,
                 power = beta, k = k, sigdf = v, sigquant = q, ntree = J,
                 ndpost = MCMCiter, nskip = burnin, keepevery = thinning)

k <- 10
L <- 20
x <- seq(min(X[, k]), max(X[, k]), length.out = L)
X.test <- X
X.test[,k] <- x[1]
for(j in 2:L){
  X.k <- X
  X.k[,k] <- x[j]
  X.test <- rbind(X.test, X.k)
} 
pred <- predict(BARTfit, X.test)
partial <- matrix(nrow = MCMCiter, ncol = L)
for(j in 1:L) {
  h <- (j - 1) * N + 1:N
  partial[, j] <- apply(pred[, h], 1, mean)
}
# Assume: partial is a 1000 x 20 matrix (1000 posterior draws, 20 grid points)
partial_df <- as.data.frame(partial)

# Add draw index for clarity
partial_df$draw <- 1:nrow(partial_df)

# Pivot to long format
long_df <- pivot_longer(partial_df, cols = -draw, names_to = "x", values_to = "value")

# Convert x to numeric (column names are likely "V1", "V2", ..., or "1", "2", ...)
long_df$x <- as.numeric(gsub("V", "", long_df$x))

summary_df <- long_df %>%
  group_by(x) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  )

ggplot(summary_df, aes(x = x, y = mean)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.4) +
  labs(x = "Grid point", y = "Partial effect", title = "Partial Dependence Function with 95% Intervals") +
  theme_minimal()
