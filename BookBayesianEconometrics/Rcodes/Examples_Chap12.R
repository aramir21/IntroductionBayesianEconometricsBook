######### Bayesian LASSO with hierarchical lambda #########
rm(list = ls()); set.seed(10101)
library(bayesreg)

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

### Using bayesreg ###

# Fit the model
fit <- bayesreg::bayesreg(y ~ X, data = df, model = "gaussian", prior = "lasso", 
                n.samples = 5000, burnin = 1000)

# Check summary
summary(fit)

# Extract posterior means of beta
beta_post_mean <- rowMeans(fit$beta)

# Compare true vs estimated
plot(beta_true, beta_post_mean, pch = 19, col = "steelblue",
     xlab = "True beta", ylab = "Posterior mean of beta",
     main = "Bayesian LASSO Shrinkage")
abline(0, 1, col = "red", lty = 2)

######### SSVS #########
rm(list = ls()); set.seed(10101)
library(BoomSpikeSlab)
library(dplyr)
library(tibble)
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
### Using BoomSpikeSlab ###
#Scale regressors
W <- scale(X); yh <- y - mean(y)
# In case you want to modify the default prior
prior <- SpikeSlabPrior(W, yh, 
                        expected.model.size = ncol(W)/2, # expect 50 nonzero predictors
                        prior.df = .01, # weaker prior than the default
                        prior.information.weight = .01,
                        diagonal.shrinkage = 0) # shrink to zero
niter <- 5000
######Estimate model########
# Default prior
SSBoomNew <- lm.spike(yh ~ W - 1, niter = niter)
Models <- SSBoomNew$beta != 0
PIP <- colMeans(SSBoomNew$beta != 0)
# Convert the logical matrix to a data frame and then to a tibble
df <- as.data.frame(Models); df_tbl <- as_tibble(df)
# Count identical rows
row_counts <- df_tbl %>% count(across(everything()), name = "frequency") %>% arrange(desc(frequency))
# Ensure your vector and matrix are logical
trueModel <- c(rep(1, 10), rep(0, 90)) == 1  # convert to logical if needed
# Assume your matrix is named 'mat'
matching_rows <- apply(row_counts[,-101], 1, function(row) all(row == trueModel))
# Get indices of best model
which(row_counts[1,] == 1)
# Frequency best model
row_counts[1,101]
# Coefficients
SummarySS <- summary(coda::mcmc(SSBoomNew$beta))
# Extract posterior means of beta
beta_post_mean <- SummarySS$statistics[, 1]
# Compare true vs estimated
plot(beta_true, beta_post_mean, pch = 19, col = "steelblue",
     xlab = "True beta", ylab = "Posterior mean of beta",
     main = "SSVS Shrinkage")
abline(0, 1, col = "red", lty = 2)

####### BART #########
rm(list = ls()); set.seed(10101)
library(BART)
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
# Train and test
c <- 0.8
Ntrain <- floor(c * N)
Ntest <- N - Ntrain
X.train <- X[1:Ntrain, ]
y.train <- y[1:Ntrain]
X.test <- X[(Ntrain+1):N, ]
y.test <- y[(Ntrain+1):N]
# Hyperparameters
alpha <- 0.95; beta <- 2; k <- 2
v <- 3; q <- 0.9; J <- 200
# MCMC parameters
MCMCiter <- 1000; burnin <- 100; thinning <- 1
# Estimate BART
BARTfit <- wbart(x.train = X.train, y.train = y.train, x.test = X.test, base = alpha,
                 power = beta, k = k, sigdf = v, sigquant = q, ntree = J,
                 ndpost = MCMCiter, nskip = burnin, keepevery = thinning)

# Trace plot sigma
library(ggplot2)
keep <- seq(burnin + 1, MCMCiter + burnin, thinning)
df_sigma <- data.frame(iteration = 1:length(keep), sigma = BARTfit$sigma[keep])
ggplot(df_sigma, aes(x = iteration, y = sigma)) +
  geom_line(color = "steelblue") +
  labs(title = "Trace Plot of Sigma", x = "Iteration", y = expression(sigma)) +
  theme_minimal()

# Prediction plot training
library(dplyr)
train_preds <- data.frame(
  y_true = y.train,
  mean = apply(BARTfit$yhat.train, 2, mean),
  lower = apply(BARTfit$yhat.train, 2, quantile, 0.025),
  upper = apply(BARTfit$yhat.train, 2, quantile, 0.975)
) %>% arrange(y_true) %>% mutate(index = row_number())
ggplot(train_preds, aes(x = index)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.4, show.legend = TRUE) +
  geom_line(aes(y = mean, color = "Predicted Mean")) +
  geom_line(aes(y = y_true, color = "True y")) +
  scale_color_manual(name = "Line", values = c("Predicted Mean" = "blue", "True y" = "black")) +
  scale_fill_manual(name = "Interval", values = c("95% Interval" = "lightblue")) +
  labs(title = "Training Data: Ordered Predictions with 95% Intervals",
       x = "Ordered Index", y = "y") +
  theme_minimal()

# Prediction plot test
test_preds <- data.frame(
  y_true = y.test,
  mean = apply(BARTfit$yhat.test, 2, mean),
  lower = apply(BARTfit$yhat.test, 2, quantile, 0.025),
  upper = apply(BARTfit$yhat.test, 2, quantile, 0.975)
) %>% arrange(y_true) %>% mutate(index = row_number())

ggplot(test_preds, aes(x = index)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Interval"), alpha = 0.4, show.legend = TRUE) +
  geom_line(aes(y = mean, color = "Predicted Mean")) +
  geom_line(aes(y = y_true, color = "True y")) +
  scale_color_manual(name = "Line", values = c("Predicted Mean" = "blue", "True y" = "black")) +
  scale_fill_manual(name = "Interval", values = c("95% Interval" = "lightblue")) +
  labs(title = "Test Data: Ordered Predictions with 95% Intervals",
       x = "Ordered Index", y = "y") +
  theme_minimal()

# Relevant regressors
Js <- c(10, 20, 50, 100, 200)
VarImportance <- matrix(0, length(Js), K)
l <- 1
for (j in Js){
  BARTfit <- wbart(x.train = X.train, y.train = y.train, x.test = X.test, base = alpha,
                   power = beta, k = k, sigdf = v, sigquant = q, ntree = j,
                   ndpost = MCMCiter, nskip = burnin, keepevery = thinning)
  VarImportance[l, ] <- BARTfit[["varcount.mean"]]/j
  l <- l + 1
}

library(tidyr)
# Assign row and column names for clarity
rownames(VarImportance) <- c("10", "20", "50", "100", "200")
colnames(VarImportance) <- as.character(1:10)

# Convert to long format for ggplot2
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
  scale_x_continuous(breaks = 1:10) +
  labs(x = "variable", y = "percent used", color = "#trees", linetype = "#trees") +
  theme_minimal()

####### Simulation: Gaussian process ########
rm(list = ls())
set.seed(10101)

library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)

# Simulation setup
n <- 100
x <- seq(0, 1, length.out = n)
sigma_f <- 1
l <- 0.2
sigma_n <- 1e-8

# Squared Exponential Kernel function
SE_kernel <- function(x1, x2, sigma_f, l) {
  outer(x1, x2, function(a, b) sigma_f^2 * exp(-0.5 * (a - b)^2 / l^2))
}

K <- SE_kernel(x, x, sigma_f, l) + diag(sigma_n, n)
samples <- mvrnorm(n = 5, mu = rep(0, n), Sigma = K)

# Transpose and rename columns to f1, f2, ..., f5
samples_t <- t(samples)
colnames(samples_t) <- paste0("f", 1:5)

# Convert to tidy data frame
df <- data.frame(x = x, samples_t) |>
  pivot_longer(cols = -x, names_to = "draw", values_to = "value")

# Plot
ggplot(df, aes(x = x, y = value, color = draw)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Simulated Gaussian Process Draws",
    x = "x", y = "f(x)", color = "Function"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

####### Given 4 observations, and one realization of a GP ###### 
rm(list = ls())
set.seed(10101)
library(ggplot2)
library(MASS)

# Define the squared exponential kernel
SE_kernel <- function(x1, x2, sigma_f, l) {
  outer(x1, x2, function(a, b) sigma_f^2 * exp(-0.5 * (a - b)^2 / l^2))
}

# Define the input space and observed points
x_star <- seq(0, 1, length.out = 200)
x0 <- c(0.1, 0.2, 0.5, 0.9)
y0 <- sin(2 * pi * x0)

# Hyperparameters
sigma_f <- 1
l <- 0.2
sigma_n <- 1e-8  # Jitter term for stability

# Compute covariance matrices
K_x0x0 <- SE_kernel(x0, x0, sigma_f, l) + diag(sigma_n, length(x0))
K_xstarx0 <- SE_kernel(x_star, x0, sigma_f, l)
K_xstarxstar <- SE_kernel(x_star, x_star, sigma_f, l) + diag(sigma_n, length(x_star))

# Compute posterior mean and covariance
K_inv <- solve(K_x0x0)
posterior_mean <- K_xstarx0 %*% K_inv %*% y0
posterior_cov <- K_xstarxstar - K_xstarx0 %*% K_inv %*% t(K_xstarx0)

# Sample from the posterior
sample_draw <- sin(2 * pi * x_star) 
# Compute 95% intervals
posterior_sd <- sqrt(diag(posterior_cov))
lower <- posterior_mean - 1.96 * posterior_sd
upper <- posterior_mean + 1.96 * posterior_sd

# Data frame for plotting
df <- data.frame(
  x = x_star,
  mean = posterior_mean,
  lower = lower,
  upper = upper,
  sample = sample_draw
)

obs <- data.frame(x = x0, y = y0)

# Plot
ggplot(df, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1.2) +
  geom_line(aes(y = sample), color = "darkgreen", linewidth = 1, linetype = "dashed") +
  geom_point(data = obs, aes(x = x, y = y), color = "red", size = 3) +
  labs(
    title = "Gaussian Process with Conditioning Points",
    x = "x", y = "f(x)",
    caption = "Blue: Posterior mean | Light blue: 95% interval | Dashed green: Population | Red: Observed points"
  ) +
  theme_minimal(base_size = 14)

####### GP with two regressors ###### 
rm(list = ls()); set.seed(10101)
# Load required packages
library(DiceKriging); library(rgl)

# Simulate training data
n_train <- 100
x1 <- runif(n_train)
x2 <- runif(n_train)
X_train <- data.frame(x1 = x1, x2 = x2)

# True function without noise
f_train <- sin(2 * pi * X_train$x1) + cos(2 * pi * X_train$x2) + sin(X_train$x1 * X_train$x2)

# Fit Gaussian Process
fit_km <- km(design = X_train, response = f_train, covtype = "gauss", nugget = 1e-10)

# Prediction grid
grid_points <- 30
x1_seq <- seq(0, 1, length.out = grid_points)
x2_seq <- seq(0, 1, length.out = grid_points)
grid <- expand.grid(x1 = x1_seq, x2 = x2_seq)

# Predict GP surface
pred <- predict(fit_km, newdata = grid, type = "UK")
z_pred <- matrix(pred$mean, nrow = grid_points, ncol = grid_points)

# Plot
persp3d(x = x1_seq, y = x2_seq, z = z_pred,
        col = "lightblue", alpha = 0.7,
        xlab = "x1", ylab = "x2", zlab = "GP Mean")
points3d(x = X_train$x1, y = X_train$x2, z = f_train, col = "red", size = 8)

fit_km@covariance@range.val # length-scale
fit_km@covariance@sd2 # Signal variance

####### Gaussian Process #########
rm(list = ls()); set.seed(10101)
library(DiceKriging); library(ggplot2)
# Simulated data
n <- 20
x <- matrix(runif(n), ncol = 1)
y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)

# Fixed values
fixed_theta <- 0.2     # length-scale
fixed_sigma2 <- 1      # process variance
fixed_nugget <- 0.01   # noise variance

# Fit model with fixed hyperparameters
fit_km <- km(
  design = data.frame(x = x),
  response = y,
  covtype = "gauss",
  nugget = fixed_nugget,
  nugget.estim = FALSE,
  coef.trend = 0,               # intercept (optional)
  coef.cov = fixed_theta,       # fixed length-scale
  coef.var = fixed_sigma2,      # fixed process variance
  estim.method = "none"         # disables hyperparameter estimation
)

# Check that parameters were fixed
fit_km@covariance@range.val  # Should return 0.2
fit_km@covariance@sd2        # Should return 1
fit_km@noise.var             # noise variance (nugget)


# Prediction grid
x_grid <- seq(0, 1, length.out = 200)
X_grid <- data.frame(x = x_grid)

# Predict mean and standard error
pred <- predict(fit_km, newdata = X_grid, type = "UK", se.compute = TRUE)

# Construct data frame for plotting
pred_df <- data.frame(
  x = x_grid,
  mean = pred$mean,
  lower = pred$mean - 1.96 * pred$sd,
  upper = pred$mean + 1.96 * pred$sd
)

# Plot
ggplot() +
  geom_ribbon(data = pred_df, aes(x = x, ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(data = pred_df, aes(x = x, y = mean), color = "blue") +
  geom_point(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "black", size = 1) +
  labs(title = "Gaussian Process Regression with DiceKriging", x = "x", y = "y") +
  theme_minimal()

fit_km@covariance@range.val  # lengthscale 
fit_km@covariance@sd2        # process variance
fit_km@noise.var             # noise variance

####### Stochastic Gradient Langevin Dynamic #########
rm(list = ls()); set.seed(10101)
library(mvtnorm)
library(MCMCpack)
library(ggplot2)
library(dplyr)

#--- Generate correlated covariates
genCovMat <- function(K, rho = 0.4) {
  Sigma0 <- matrix(1, K, K)
  for (i in 2:K) {
    for (j in 1:(i - 1)) {
      Sigma0[i, j] <- runif(1, -rho, rho)^(i - j)
    }
  }
  Sigma0 <- Sigma0 * t(Sigma0)
  diag(Sigma0) <- 1
  return(Sigma0)
}

#--- Simulate logistic regression data
simulate_logit_data <- function(K, N, beta_true) {
  Sigma0 <- genCovMat(K)
  X <- rmvnorm(N, mean = rep(0, K), sigma = Sigma0)
  linpred <- X %*% beta_true
  p <- 1 / (1 + exp(-linpred))
  y <- rbinom(N, 1, p)
  list(y = y, X = X)
}

#--- One SGLD step
SGLD_step <- function(beta, y, X, stepsize, batch_size, prior_var = 10) {
  N <- nrow(X)
  K <- length(beta)
  ids <- sample(1:N, size = batch_size, replace = FALSE)
  
  grad <- rep(0, K)
  for (i in ids) {
    xi <- X[i, ]
    eta <- sum(xi * beta)
    pi <- 1 / (1 + exp(-eta))
    grad_i <- -(y[i] - pi) * xi
    grad <- grad + grad_i
  }
  
  grad <- grad / batch_size * N
  grad <- grad + beta / prior_var  # gradient of log-prior
  
  noise <- rnorm(K, 0, sqrt(stepsize))
  beta_new <- beta - 0.5 * stepsize * grad + noise
  return(beta_new)
}

#--- SGLD algorithm
run_SGLD <- function(y, X, stepsize, batch_prop, n_iter, burnin, beta_init = NULL) {
  N <- nrow(X)
  K <- ncol(X)
  batch_size <- round(batch_prop * N)
  
  beta_mat <- matrix(0, n_iter + burnin, K)
  beta_mat[1, ] <- if (is.null(beta_init)) rep(0, K) else beta_init
  
  for (s in 2:(n_iter + burnin)) {
    beta_mat[s, ] <- SGLD_step(beta_mat[s - 1, ], y, X, stepsize, batch_size)
  }
  beta_mat[(burnin + 1):(n_iter + burnin), ]
}

#--- Parameters
K <- 10
N <- 100000
beta_true <- rep(0.5, K)
batch_prop <- 0.01
n_iter <- 2000
burnin <- 500
stepsize <- 1e-4
k_target <- 5  # beta5

#--- Simulate data
sim_data <- simulate_logit_data(K, N, beta_true)
y <- sim_data$y
X <- sim_data$X

#--- Run SGLD
posterior_sgld <- run_SGLD(y, X, stepsize, batch_prop, n_iter, burnin)

#--- Run MCMCpack logit
df <- as.data.frame(X)
colnames(df) <- paste0("X", 1:K)
df$y <- y
formula <- as.formula(paste("y ~", paste(colnames(df)[1:K], collapse = " + "), "-1"))
posterior_mh <- MCMClogit(formula, data = df, b0 = 0, B0 = 0.1,
                          burnin = burnin, mcmc = n_iter)

#--- Compare densities for beta5
df_plot <- data.frame(
  value = c(posterior_sgld[, k_target], posterior_mh[, k_target]),
  method = rep(c("SGLD", "MCMC"), each = n_iter)
)

ggplot(df_plot, aes(x = value, fill = method, color = method)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = beta_true[k_target], linetype = "dashed", color = "black") +
  labs(title = expression(paste("Posterior density of ", beta[5])),
       x = expression(beta[5]),
       y = "Density") +
  theme_minimal()
