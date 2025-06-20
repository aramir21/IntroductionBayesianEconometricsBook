######### Bayesian LASSO with hierarchical lambda
rm(list = ls()); set.seed(10101)
library(monomvn)

# Parameters
n <- 500  # sample size
p <- 100  # number of predictors
s <- 10   # number of non-zero coefficients

# Generate design matrix
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# True beta: first s coefficients are non-zero, rest are zero
beta_true <- c(runif(s, -3, 3), rep(0, p - s))

# Generate response with some noise
sigma <- 1
y <- X %*% beta_true + rnorm(n, sd = sigma)
df <- data.frame(X,y)

### Using monomvn ###
library(monomvn)

# Fit model using Bayesian LASSO
fit <- monomvn::blasso(X, y, T = 5000, verb = 0)

# burn-in removal
burnin <- 1000
beta_samples <- fit$beta[(burnin + 1):5000, ]  # exclude burn-in
intercept_samples <- fit$mu[(burnin + 1):5000] # intercept posterior draws

# Posterior mean of coefficients
beta_post_mean <- colMeans(beta_samples)

# Plot true vs estimated betas
plot(beta_true, beta_post_mean, pch = 19, col = "steelblue",
     xlab = "True beta", ylab = "Posterior mean of beta",
     main = "Bayesian LASSO (blasso): Shrinkage Effect")
abline(0, 1, col = "red", lty = 2)

# Intercept summary
mean(intercept_samples)
quantile(intercept_samples, c(0.025, 0.975))

# Local shrinkage parameters
colMeans(is.na(fit$tau2i))  # proportion of NA values per coefficient

