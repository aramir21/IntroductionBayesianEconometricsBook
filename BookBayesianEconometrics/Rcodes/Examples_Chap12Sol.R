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

######### Application BART: Conflict #########
rm(list = ls()); set.seed(10101)
library(BART); library(tidyr)
library(pROC); library(caret)

# Load and prepare data
Data <- read.csv("https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/Conflict.csv", sep = ",", header = TRUE, quote = "")
X <- as.matrix(Data[, -c(1, 2)])
y <- unlist(Data[, 2])

niter <- 1000; burnin <- 200; J <- 200

# Create 5-fold cross-validation indices
folds <- createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)

thresholds <- seq(0.01, 0.99, by = 0.01)
scores <- matrix(0, nrow = 5, ncol = length(thresholds))

for (i in 1:5) {
  test_idx <- folds[[i]]
  train_idx <- setdiff(1:length(y), test_idx)
  
  X.train <- X[train_idx, ]
  y.train <- y[train_idx]
  X.test <- X[test_idx, ]
  y.test <- y[test_idx]
  
  pf <- pbart(x.train = X.train, y = y.train, x.test = X.train, ntree = J, ndpost = niter, nskip = burnin)
  prob_train <- pf$prob.test.mean
  
  for (j in 1:length(thresholds)) {
    thresh <- thresholds[j]
    preds <- ifelse(prob_train > thresh, 1, 0)
    TP <- sum(preds == 1 & y.train == 1)
    TN <- sum(preds == 0 & y.train == 0)
    FP <- sum(preds == 1 & y.train == 0)
    FN <- sum(preds == 0 & y.train == 1)
    TPR <- TP / (TP + FN + 1e-8)
    TNR <- TN / (TN + FP + 1e-8)
    scores[i, j] <- TPR + TNR
  }
}

mean_scores <- colMeans(scores)
best_thresh <- thresholds[which.max(mean_scores)]
cat("Optimal threshold by CV:", best_thresh, "\n")

# Train/test split to evaluate final performance
N <- nrow(X); K <- ncol(X)
N.train <- floor(0.8 * N)
X.train <- X[1:N.train, ]
y.train <- y[1:N.train]
X.test <- X[(N.train + 1):N, ]
y.test <- y[(N.train + 1):N]

pf_final <- pbart(x.train = X.train, y = y.train, x.test = X.test, ntree = 200, ndpost = 1000, nskip = 200)
prob_test <- pf_final$prob.test.mean
y_pred_test <- ifelse(prob_test > best_thresh, 1, 0)

conf_mat <- table(Predicted = y_pred_test, Actual = y.test)
accuracy <- mean(y_pred_test == y.test)

print(conf_mat)
cat("Test set accuracy:", round(accuracy, 3), "\n")

# Relevant regressors
Js <- c(10, 20, 50, 100, 200)
VarImportance <- matrix(0, length(Js), K)
l <- 1
for (j in Js){
  BARTfit <- BART::pbart(x.train = X.train, y = y.train, x.test = X.test, ntree = j, ndpost = niter, nskip = burnin)
  VarImportance[l, ] <- BARTfit[["varcount.mean"]]/j
  l <- l + 1
}

# Assign row and column names for clarity
rownames(VarImportance) <- as.character(Js)
colnames(VarImportance) <- as.character(1:K)

# Convert to long format
importance_df <- as.data.frame(VarImportance) %>%
  mutate(trees = rownames(.)) %>%
  pivot_longer(
    cols = -trees,
    names_to = "variable",
    values_to = "percent_used"
  )

# Make variable numeric for plotting
importance_df$variable <- as.numeric(importance_df$variable)

importance_df$trees <- factor(importance_df$trees, levels = c("10", "20", "50", "100", "200"))

ggplot(importance_df, aes(x = variable, y = percent_used, color = trees, linetype = trees)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("10" = "red", "20" = "green", "50" = "blue", "100" = "cyan", "200" = "magenta")) +
  scale_x_continuous(breaks = 1:K) +
  labs(x = "variable", y = "percent used", color = "#trees", linetype = "#trees") +
  theme_minimal()

data.frame(Nanes = colnames(X), Frequency = VarImportance[1,])

######## Gaussian Process #########
rm(list = ls()); set.seed(10101)
library(fields); library(DiceKriging); library(ggplot2)
# Simulate training data
n_train <- 500
x1 <- runif(n_train)
x2 <- runif(n_train)
X_train <- data.frame(x1 = x1, x2 = x2)
sig <- 0.1
u <- rnorm(n_train, mean = 0, sd = sig)
# True function without noise
f_train <- sin(2 * pi * X_train$x1) + cos(2 * pi * X_train$x2) + sin(X_train$x1 * X_train$x2)
y_train <- scale(f_train + u)

# Empirical Bayes: Optimize the log marginal likelihood
log_marginal_likelihood <- function(par, X, y) {
  sigma2_f <- par[1]
  l2 <- par[2]
  sigma2_n <- par[3]
  n <- nrow(X)
  
  # Ensure parameters are positive
  if (sigma2_f <= 0 | l2 <= 0 | sigma2_n <= 0) return(Inf)
  
  # Squared exponential kernel
  dists <- rdist(X)
  K <- sigma2_f * exp(-0.5 * (dists^2) / l2)
  
  # Add noise variance + jitter for numerical stability
  K <- K + (sigma2_n + 1e-8) * diag(n)
  
  # Cholesky decomposition
  L <- tryCatch(chol(K), error = function(e) return(NULL))
  if (is.null(L)) return(Inf)
  
  alpha <- backsolve(t(L), forwardsolve(L, y)) #--> K^{-1}y
  log_det_K <- 2 * sum(log(diag(L)))
  lml <- -0.5 * t(y) %*% alpha - 0.5 * log_det_K - 0.5 * n * log(2 * pi)
  
  return(-as.numeric(lml))  # Return negative log-marginal likelihood
}
par0 <- rep(1, 3)
# log_marginal_likelihood(par = par0, X = X_train, y = y_train)
# log_marginal_likelihood(par = c(11.6, 0.52^2, 0.1^2), X = X_train, y = y_train)
ResOpt <- optim(par = par0, fn = log_marginal_likelihood, method = "L-BFGS-B",
                lower = c(1e-5, 1e-5), X = X_train, y = y_train)
# ResOpt$par
# ResOpt$value
theta_hat <- ResOpt$par
sigma2_f_hat <- theta_hat[1]
l2_hat <- theta_hat[2]
sigma2_n_hat <- theta_hat[3]

# Grid for prediction 
grid_points <- 20
x1_seq <- seq(0, 1, length.out = grid_points)
x2_seq <- seq(0, 1, length.out = grid_points)
X_new <- expand.grid(x1 = x1_seq, x2 = x2_seq)

# Fit Gaussian Process
fit_km <- km(design = X_train, response = y_train, covtype = "gauss",
             noise.var = rep(sigma2_n_hat, n_train))

# Predict GP surface
pred <- predict(fit_km, newdata = X_new, type = "UK")
z_pred <- matrix(pred$mean, nrow = grid_points, ncol = grid_points)
# Plot
persp3d(x = x1_seq, y = x2_seq, z = z_pred,
        col = "lightblue", alpha = 0.7,
        xlab = "x1", ylab = "x2", zlab = "GP Mean")
points3d(x = X_train$x1, y = X_train$x2, z = y_train, col = "red", size = 8)

# fit_km@covariance@range.val # length-scale
# fit_km@covariance@sd2 # Signal variance
# log_marginal_likelihood(par = c(fit_km@covariance@sd2, fit_km@covariance@range.val^2, sigma2_n_hat), X = X_train, y = y_train)


# # Covariance matrices
# K <- sigma2_f_hat * exp(-0.5 * rdist(X_train)^2 / l2_hat) +
#   (sigma2_n_hat + 1e-4) * diag(n_train)
# 
# K_s <- sigma2_f_hat * exp(-0.5 * rdist(X_new, X_train)^2 / l2_hat)
# K_ss <- sigma2_f_hat * exp(-0.5 * rdist(X_new)^2 / l2_hat)
# 
# # Prediction new points
# L <- chol(K)
# alpha <- backsolve(t(L), forwardsolve(L, y_train)) # --> K^{-1}y
# 
# # Posterior mean
# mu_post <- K_s %*% alpha
# z_pred <- matrix(mu_post, nrow = grid_points, ncol = grid_points)
# 
# # Posterior covariance
# v <- forwardsolve(L, t(K_s)) # --> L^{-1}K_s'
# cov_post <- K_ss - t(v) %*% v
# sd_post <- sqrt(pmax(diag(cov_post), 0))  # ensure no negative values
# 
# # Plot
# persp3d(x = x1_seq, y = x2_seq, z = z_pred,
#         col = "lightblue", alpha = 0.7,
#         xlab = "x1", ylab = "x2", zlab = "GP Mean")
# points3d(x = X_train$x1, y = X_train$x2, z = y_train, col = "red", size = 8)

####### Stochastic Gradient Langevin Dynamic #########
rm(list = ls()); set.seed(10101)
