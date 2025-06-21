######### Bayesian LASSO with hierarchical lambda
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
