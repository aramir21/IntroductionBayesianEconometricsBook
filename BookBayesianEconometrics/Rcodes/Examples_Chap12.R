######### Bayesian LASSO with hierarchical lambda
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

