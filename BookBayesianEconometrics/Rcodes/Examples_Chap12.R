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
