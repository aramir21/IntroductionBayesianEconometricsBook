########################## Simulation exercise: Gaussian mixture with 2 components no regressors########################## 
rm(list = ls())
set.seed(010101)
library(ggplot2)

# Simulate data from a 2-component mixture model
n <- 500
z <- rbinom(n, 1, 0.75)  # Latent class indicator
y <- ifelse(z == 0, rnorm(n, 0.5, 1), rnorm(n, 2.5, 1))
data <- data.frame(y)

# Plot
ggplot(data, aes(x = y)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(title = "Density Plot", x = "y", y = "Density") +
  theme_minimal()

# Hyperparameters
mu0 <- 0
sig2mu0 <- 10
H <- 2
a0h <- rep(1/H, H)
# MCMC parameters
mcmc <- 1000
burnin <- 500
tot <- mcmc + burnin
thin <- 2
# Gibbs sampling functions
Postmu <- function(yh){
  Nh <- length(yh)
  sig2mu <- (1/sig2mu0 + Nh)^(-1)
  mun <- sig2mu*(mu0/sig2mu0 + sum(yh))
  mu <- rnorm(1, mun, sig2mu^0.5)
  return(mu)
}

PostPsi <- matrix(NA, tot, n)
PostMu <- matrix(NA, tot, H)
PostLambda <- rep(NA, tot)
Id1 <- which(y <= 1) # 1 is from inspection of the density plot of y 
Id2 <- which(y > 1)
N1 <- length(Id1); N2 <- length(Id2)
Lambda <- c(N1/n, N2/n)
MU <- c(mean(y[Id1]), mean(y[Id2]))
Psi <- rep(NA, n)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  for(i in 1:n){
    lambdai <- NULL
    for(h in 1:H){
      lambdaih <- Lambda[h]*dnorm(y[i], MU[h], 1)
      lambdai <- c(lambdai, lambdaih)
    }
    Psi[i] <- sample(1:H, 1, prob = lambdai)
  }
  PostPsi[s, ] <- Psi
  for(h in 1:H){
    idh <- which(Psi == h)
    MU[h] <- Postmu(yh = y[idh])
  }
  PostMu[s,] <- MU 
  Lambda <- MCMCpack::rdirichlet(1, a0h + table(Psi))
  PostLambda[s] <- Lambda[1]
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq(burnin, tot, thin)
PosteriorMUs <- coda::mcmc(PostMu[keep,])
summary(PosteriorMUs)
plot(PosteriorMUs)
dfMU <- data.frame(mu1 = PostMu[keep,1], mu2 = PostMu[keep,2])
# Plot
require(latex2exp)
ggplot(dfMU) +
  geom_density(aes(x = mu1, color = "mu1"), linewidth = 1) +  # First density plot
  geom_density(aes(x = mu2, color = "mu2"), linewidth = 1) +  # Second density plot
  labs(title = "Density Plot", x = TeX("$\\mu$"), y = "Density", color = "Variable") +
  theme_minimal() +
  scale_color_manual(values = c("mu1" = "blue", "mu2" = "red"))  # Custom colors

PosteriorLAMBDA <- coda::mcmc(PostLambda[keep])
summary(PosteriorLAMBDA)
plot(PosteriorLAMBDA)

##### Change mean #####
set.seed(010101)
library(ggplot2)

# Simulate data from a 2-component mixture model
n <- 500
z <- rbinom(n, 1, 0.75)  # Latent class indicator
y <- ifelse(z == 0, rnorm(n, 0.5, 1), rnorm(n, 1, 1))
data <- data.frame(y)

# Plot
ggplot(data, aes(x = y)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(title = "Density Plot", x = "y", y = "Density") +
  theme_minimal()

# Hyperparameters
mu0 <- 0
sig2mu0 <- 10
H <- 2
a0h <- rep(1/H, H)
# MCMC parameters
mcmc <- 1000
burnin <- 500
tot <- mcmc + burnin
thin <- 2
# Gibbs sampling functions
Postmu <- function(yh){
  Nh <- length(yh)
  sig2mu <- (1/sig2mu0 + Nh)^(-1)
  mun <- sig2mu*(mu0/sig2mu0 + sum(yh))
  mu <- rnorm(1, mun, sig2mu^0.5)
  return(mu)
}

PostPsi <- matrix(NA, tot, n)
PostMu <- matrix(NA, tot, H)
PostLambda <- rep(NA, tot)
Id1 <- which(y <= 1) # 1 is from inspection of the density plot of y 
Id2 <- which(y > 1)
N1 <- length(Id1); N2 <- length(Id2)
Lambda <- c(N1/n, N2/n)
MU <- c(mean(y[Id1]), mean(y[Id2]))
Psi <- rep(NA, n)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  for(i in 1:n){
    lambdai <- NULL
    for(h in 1:H){
      lambdaih <- Lambda[h]*dnorm(y[i], MU[h], 1)
      lambdai <- c(lambdai, lambdaih)
    }
    Psi[i] <- sample(1:H, 1, prob = lambdai)
  }
  PostPsi[s, ] <- Psi
  for(h in 1:H){
    idh <- which(Psi == h)
    MU[h] <- Postmu(yh = y[idh])
  }
  PostMu[s,] <- MU 
  Lambda <- MCMCpack::rdirichlet(1, a0h + table(Psi))
  PostLambda[s] <- Lambda[1]
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq(burnin, tot, thin)
PosteriorMUs <- coda::mcmc(PostMu[keep,])
summary(PosteriorMUs)
plot(PosteriorMUs)
dfMU <- data.frame(mu1 = PostMu[keep,1], mu2 = PostMu[keep,2])
# Plot
require(latex2exp)
ggplot(dfMU) +
  geom_density(aes(x = mu1, color = "mu1"), linewidth = 1) +  # First density plot
  geom_density(aes(x = mu2, color = "mu2"), linewidth = 1) +  # Second density plot
  labs(title = "Density Plot", x = TeX("$\\mu$"), y = "Density", color = "Variable") +
  theme_minimal() +
  scale_color_manual(values = c("mu1" = "blue", "mu2" = "red"))  # Custom colors

PosteriorLAMBDA <- coda::mcmc(PostLambda[keep])
summary(PosteriorLAMBDA)
plot(PosteriorLAMBDA)

###### Permutations ######
rm(list = ls())
set.seed(010101)
library(ggplot2)

# Simulate data from a 2-component mixture model
n <- 500
z <- rbinom(n, 1, 0.75)  # Latent class indicator
y <- ifelse(z == 0, rnorm(n, 0.5, 1), rnorm(n, 1, 1))
data <- data.frame(y)

# Plot
ggplot(data, aes(x = y)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(title = "Density Plot", x = "y", y = "Density") +
  theme_minimal()

# Hyperparameters
mu0 <- 0
sig2mu0 <- 10
H <- 2
a0h <- rep(1/H, H)
# MCMC parameters
mcmc <- 2000
burnin <- 500
tot <- mcmc + burnin
thin <- 2
# Gibbs sampling functions
Postmu <- function(yh){
  Nh <- length(yh)
  sig2mu <- (1/sig2mu0 + Nh)^(-1)
  mun <- sig2mu*(mu0/sig2mu0 + sum(yh))
  mu <- rnorm(1, mun, sig2mu^0.5)
  return(mu)
}

PostPsi <- matrix(NA, tot, n)
PostMu <- matrix(NA, tot, H)
PostLambda <- rep(NA, tot)
Id1 <- which(y <= 1) # 1 is from inspection of the density plot of y 
Id2 <- which(y > 1)
N1 <- length(Id1); N2 <- length(Id2)
Lambda <- c(N1/n, N2/n)
MU <- c(mean(y[Id1]), mean(y[Id2]))
Psi <- rep(NA, n)
per1 <- c(1,2)
per2 <- c(2,1)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  for(i in 1:n){
    lambdai <- NULL
    for(h in 1:H){
      lambdaih <- Lambda[h]*dnorm(y[i], MU[h], 1)
      lambdai <- c(lambdai, lambdaih)
    }
    Psi[i] <- sample(1:H, 1, prob = lambdai)
  }
  for(h in 1:H){
    idh <- which(Psi == h)
    MU[h] <- Postmu(yh = y[idh])
  }
  Lambda <- MCMCpack::rdirichlet(1, a0h + table(Psi))
  # Permutations
  labels <- sample(1:2, 1, prob = c(0.5, 0.5))
  if(labels == 2){
    Lambda <- Lambda[per2]
    MU <- MU[per2]
    for(i in 1:n){
      if(Psi[i] == 1){
        Psi[i] <- 2
      }else{
        Psi[i] <- 1
      }
    }
  }
  PostPsi[s, ] <- Psi
  PostMu[s,] <- MU 
  PostLambda[s] <- Lambda[1]
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq(burnin, tot, thin)
PosteriorMUs <- coda::mcmc(PostMu[keep,])
summary(PosteriorMUs)
plot(PosteriorMUs)
dfMU <- data.frame(mu1 = PostMu[keep,1], mu2 = PostMu[keep,2])
# Plot
require(latex2exp)
ggplot(dfMU) +
  geom_density(aes(x = mu1, color = "mu1"), linewidth = 1) +  # First density plot
  geom_density(aes(x = mu2, color = "mu2"), linewidth = 1) +  # Second density plot
  labs(title = "Density Plot", x = TeX("$\\mu$"), y = "Density", color = "Variable") +
  theme_minimal() +
  scale_color_manual(values = c("mu1" = "blue", "mu2" = "red"))  # Custom colors

PosteriorLAMBDA <- coda::mcmc(PostLambda[keep])
summary(PosteriorLAMBDA)
plot(PosteriorLAMBDA)

########################## Simulation exercise: Gaussian mixture: 2 components ########################## 
rm(list = ls())
set.seed(010101)
library(brms)
library(ggplot2)

# Simulate data from a 2-component mixture model
n <- 1000
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
rm(list = ls())
set.seed(010101)
# Simulate data from a 2-component mixture model
n <- 1000
x <- rnorm(n)
z <- rbinom(n, 1, 0.5)  # Latent class indicator
y <- ifelse(z == 0, rnorm(n, 2 + 1.5*x, 1), rnorm(n, -1 + 0.5*x, 0.8))
data <- data.frame(y, x)
# Hyperparameters
d0 <- 0.001
a0 <- 0.001
b0 <- rep(0, 2)
B0 <- diag(2)
B0i <- solve(B0)
a01 <- 1/2
a02 <- 1/2
# MCMC parameters
mcmc <- 5000
burnin <- 1000
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

PostBetas1 <- matrix(0, tot, 2)
PostBetas2 <- matrix(0, tot, 2)
PostSigma21 <- rep(0, tot)
PostSigma22 <- rep(0, tot)
PostPsi <- matrix(0, tot, n)
PostLambda <- rep(0, tot)
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
PosteriorLambda1 <- coda::mcmc(PostLambda[keep])
summary(PosteriorLambda1)
plot(PosteriorLambda1)

########################## Simulation exercise: Semi-parametric ########################## 
rm(list = ls())
set.seed(010101)
library(ggplot2)
# Simulate data from a 2-component mixture model
n <- 500
x1 <- rnorm(n); x2 <- rnorm(n)
X <- cbind(x1,x2); B <- c(-0.5, 1.5)
u <- rt(n, 3)
y <- 1 + X%*%B + u
Reg <- lm(y ~ X)
Res <- Reg$residuals
data <- data.frame(Res)
# Plot
ggplot(data, aes(x = Res)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(title = "Density Plot", x = "Residuals", y = "Density") +
  theme_minimal()

# Hyperparameters
d0 <- 0.001
a0 <- 0.001
b0 <- rep(0, 2)
B0 <- diag(2)
B0i <- solve(B0)
mu0 <- 0
sig2mu0 <- 10
H <- 5
a0h <- rep(1/H, H)
# MCMC parameters
mcmc <- 2000
burnin <- 4000
tot <- mcmc + burnin
thin <- 2
# Gibbs sampling functions
PostSig2 <- function(Beta, muh, Xh, yh){
  Nh <- length(yh)
  an <- a0 + Nh
  dn <- d0 + t(yh - muh - Xh%*%Beta)%*%(yh - muh - Xh%*%Beta)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBeta <- function(sig2, mu, X, y, Psi){
  XtX <- matrix(0, 2, 2)
  Xty <- matrix(0, 2, 1)
  Hs <- length(mu)
  for(h in 1:Hs){
    idh <- which(Psi == h)
    if(length(idh) == 1){
      Xh <- matrix(X[idh,], 1, 2)
      XtXh <- sig2[h]^(-1)*t(Xh)%*%Xh
      yh <- y[idh]
      Xtyh <- sig2[h]^(-1)*t(Xh)%*%(yh - mu[h])
    }else{
      Xh <- X[idh,]
      XtXh <- sig2[h]^(-1)*t(Xh)%*%Xh
      yh <- y[idh]
      Xtyh <- sig2[h]^(-1)*t(Xh)%*%(yh - mu[h])
    }
    XtX <- XtX + XtXh
    Xty <- Xty + Xtyh
  }
  Bn <- solve(B0i + XtX)
  bn <- Bn%*%(B0i%*%b0 + Xty)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
Postmu <- function(sig2h, Beta, Xh, yh){
  Nh <- length(yh)
  sig2mu <- (1/sig2mu0 + Nh/sig2h)^(-1)
  mun <- sig2mu*(mu0/sig2mu0 + sum((yh - Xh%*%Beta))/sig2h)
  mu <- rnorm(1, mun, sig2mu^0.5)
  return(mu)
}

PostBetas <- matrix(0, tot, 2)
PostPsi <- matrix(0, tot, n)
PostSigma2 <- list()
PostMu <- list()
PostLambda <- list()
Resq <- quantile(Res, c(0.2, 0.4, 0.6, 0.8))
Id1 <- which(Res <= Resq[1])
Id2 <- which(Res > Resq[1] & Res <= Resq[2])
Id3 <- which(Res > Resq[2] & Res <= Resq[3])
Id4 <- which(Res > Resq[3] & Res <= Resq[4])
Id5 <- which(Res > Resq[4])
Nh <- rep(n/H, H)
Lambda <- rep(1/H, H)
MU <- c(mean(Res[Id1]), mean(Res[Id2]), mean(Res[Id3]), mean(Res[Id4]), mean(Res[Id5]))
Sig2 <- c(var(Res[Id1]), var(Res[Id2]), var(Res[Id3]), var(Res[Id4]), var(Res[Id5]))
Beta <- Reg$coefficients[2:3]
Psi <- rep(NA, n)
Hs <- length(MU)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  for(i in 1:n){
    lambdai <- NULL
    for(h in 1:Hs){
      lambdaih <- Lambda[h]*dnorm(y[i] - X[i,]%*%Beta, MU[h], Sig2[h]^0.5)
      lambdai <- c(lambdai, lambdaih)
    }
    Psi[i] <- sample(1:Hs, 1, prob = lambdai)
  }
  PostPsi[s, ] <- Psi
  Hs <- length(table(Psi))
  for(h in 1:Hs){
    idh <- which(Psi == h)
    Sig2[h] <- PostSig2(Beta = Beta, muh = MU[h], Xh = X[idh,], yh = y[idh])
    MU[h] <- Postmu(sig2h = Sig2[h], Beta = Beta, Xh = X[idh,], yh = y[idh])
  }
  PostSigma2[[s]] <- Sig2
  PostMu[[s]] <- MU 
  Beta <- PostBeta(sig2 = Sig2, mu = MU, X = X, y = y, Psi = Psi)
  PostBetas[s,] <- Beta
  Lambda <- sort(MCMCpack::rdirichlet(1, a0h[1:Hs] + table(Psi)), decreasing = TRUE)
  PostLambda[[s]] <- Lambda
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq(burnin, tot, thin)
PosteriorBetas <- coda::mcmc(PostBetas[keep,])
summary(PosteriorBetas)
plot(PosteriorBetas)

PosteriorPsi <- PostPsi[keep,]
Clusters <- sapply(1:length(keep), function(i){length(table(PosteriorPsi[i,]))})
NClus <- 2
PosteriorSIGMA <- matrix(NA, length(keep), NClus)
PosteriorMU <- matrix(NA, length(keep), NClus)
PosteriorLAMBDA <- matrix(NA, length(keep), NClus)
l <- 1
for (s in keep){
  PosteriorSIGMA[l,] <- PostSigma2[[s]][1:NClus]
  PosteriorMU[l,] <- PostMu[[s]][1:NClus]
  PosteriorLAMBDA[l,] <- PostLambda[[s]][1:NClus]
  l <- l + 1
}

summary(coda::mcmc(PosteriorSIGMA))
summary(coda::mcmc(PosteriorMU))
summary(coda::mcmc(PosteriorLAMBDA))
