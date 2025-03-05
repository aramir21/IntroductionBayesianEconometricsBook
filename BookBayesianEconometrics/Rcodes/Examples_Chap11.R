########################## Simulation exercise: Gaussian mixture: 2 components ########################## 
rm(list = ls())
set.seed(010101)
library(brms)
library(ggplot2)

# Simulate data from a 2-component mixture model
n <- 500
x <- rnorm(n)
z <- rbinom(n, 1, 0.5)  # Latent class indicator
y <- ifelse(z == 0, rnorm(n, 2 + 1.5*x, 1), rnorm(n, -1 + 0.5*x, 0.8))
data <- data.frame(y, x)

# Plot
ggplot(data, aes(x = y)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(title = "Density Plot", x = "y", y = "Density") +
  theme_minimal()

# Define priors
priors <- c(
  set_prior("normal(0, 5)", class = "Intercept", dpar = "mu1"),  # First component intercept
  set_prior("normal(0, 5)", class = "Intercept", dpar = "mu2"),  # Second component intercept
  set_prior("normal(0, 2)", class = "b", dpar = "mu1"),  # First component slope
  set_prior("normal(0, 2)", class = "b", dpar = "mu2"),  # Second component slope
  set_prior("cauchy(0, 2)", class = "sigma1", lb = 0),  # First component sigma
  set_prior("cauchy(0, 2)", class = "sigma2", lb = 0),  # Second component sigma
  set_prior("dirichlet(1, 1)", class = "theta")  # Mixing proportions
)

# Fit a 2-component Gaussian mixture regression model
fit <- brm(
  bf(y ~ 1 + x, family = mixture(gaussian, gaussian)),  # Two normal distributions
  data = data,
  prior = priors,
  chains = 4, iter = 2000, warmup = 1000, cores = 4
)
prior_summary(fit) # Summary of priors
summary(fit) # Summary of posterior draws
plot(fit) # Plots of posterior draws

########### Perform inference from scratch ###############
# Hyperparameters
d0 <- 0.001
a0 <- 0.001
b0 <- rep(0, 2)
B0 <- diag(2)
B0i <- solve(B0)
a01 <- 1/2
a02 <- 1/2
# MCMC parameters
mcmc <- 2000
burnin <- 500
tot <- mcmc + burnin
thin <- 2
# Gibbs sampling functions
PostSig2 <- function(Betah, Xh, yh){
  Nh <- length(yh)
  an <- a0 + Nh
  dn <- d0 + t(yh - Xh%*%Betah)%*%(yh - Xh%*%Betah)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBeta <- function(sig2h, Xh, yh){
  Bn <- solve(B0i + sig2h^(-1)*t(Xh)%*%Xh)
  bn <- Bn%*%(B0i%*%b0 + sig2h^(-1)*t(Xh)%*%yh)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}

PostBetas1 <- matrix(0, mcmc+burnin, 2)
PostBetas2 <- matrix(0, mcmc+burnin, 2)
PostSigma21 <- rep(0, mcmc+burnin)
PostSigma22 <- rep(0, mcmc+burnin)
PostPsi <- matrix(0, mcmc+burnin, n)
PostLambda <- rep(0, mcmc+burnin)
Id1 <- which(y<1) # 1 is from inspection of the density plot of y 
N1 <- length(Id1)
Lambda1 <- N1/n
Id2 <- which(y>=1) # 1 is from inspection of the density plot of y 
N2 <- length(Id2)
Lambda2 <- N2/n
Reg1 <- lm(y ~ x, subset = Id1)
SumReg1 <- summary(Reg1)
Beta1 <- Reg1$coefficients
sig21 <- SumReg1$sigma^2 
Reg2 <- lm(y ~ x, subset = Id2)
SumReg2 <- summary(Reg2)
Beta2 <- Reg2$coefficients
sig22 <- SumReg2$sigma^2
X <- cbind(1, x)
Psi <- rep(NA, n)
for(s in 1:tot){
  for(i in 1:n){
    lambdai1 <- Lambda1*dnorm(y[i], X[i,]%*%Beta1, sig21^0.5)
    lambdai2 <- Lambda2*dnorm(y[i], X[i,]%*%Beta2, sig22^0.5)
    Psi[i] <- sample(c(1,2), 1, prob = c(lambdai1, lambdai2))
  }
  PostPsi[s, ] <- Psi
  Id1 <- which(Psi == 1); Id2 <- which(Psi == 2)
  N1 <- length(Id1); N2 <- length(Id2)
  sig21 <- PostSig2(Betah = Beta1, Xh = X[Id1, ], yh = y[Id1])
  sig22 <- PostSig2(Betah = Beta2, Xh = X[Id2, ], yh = y[Id2])
  PostSigma21[s] <- sig21
  PostSigma22[s] <- sig22
  Beta1 <- PostBeta(sig2h = sig21, Xh = X[Id1, ], yh = y[Id1])
  Beta2 <- PostBeta(sig2h = sig22, Xh = X[Id2, ], yh = y[Id2])
  PostBetas1[s,] <- Beta1
  PostBetas2[s,] <- Beta2
  Lambda <- MCMCpack::rdirichlet(1, c(a01 + N1, a02 + N2))
  Lambda1 <- Lambda[1]; Lambda2 <- Lambda[2]
  PostLambda[s] <- Lambda1 
}
keep <- seq((burnin+1), tot, thin)
PosteriorBetas1 <- coda::mcmc(PostBetas1[keep,])
summary(PosteriorBetas1)
plot(PosteriorBetas1)
PosteriorBetas2 <- coda::mcmc(PostBetas2[keep,])
summary(PosteriorBetas2)
plot(PosteriorBetas2)
PosteriorSigma21 <- coda::mcmc(PostSigma21[keep])
summary(PosteriorSigma21)
plot(PosteriorSigma21)
PosteriorSigma22 <- coda::mcmc(PostSigma22[keep])
summary(PosteriorSigma22)
plot(PosteriorSigma22)

