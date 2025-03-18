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
  Lambda <- sort(MCMCpack::rdirichlet(1, a0h + table(Psi)), decreasing = TRUE)
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
  Lambda <- sort(MCMCpack::rdirichlet(1, a0h + table(Psi)), decreasing = TRUE)
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
rm(list = ls()); set.seed(010101); library(ggplot2)
# Simulate data from a 2-component mixture model
n <- 500
z <- rbinom(n, 1, 0.75)  # Latent class indicator
y <- ifelse(z == 0, rnorm(n, 0.5, 1), rnorm(n, 2.5, 1))
# Hyperparameters
mu0 <- 0; sig2mu0 <- 10; H <- 2; a0h <- rep(1/H, H)
# MCMC parameters
mcmc <- 2000; burnin <- 500
tot <- mcmc + burnin; thin <- 2
# Gibbs sampling functions
Postmu <- function(yh){
  Nh <- length(yh)
  sig2mu <- (1/sig2mu0 + Nh)^(-1)
  mun <- sig2mu*(mu0/sig2mu0 + sum(yh))
  mu <- rnorm(1, mun, sig2mu^0.5)
  return(mu)
}
PostPsi <- matrix(NA, tot, n); PostMu <- matrix(NA, tot, H)
PostLambda <- rep(NA, tot)
Id1 <- which(y <= 1); Id2 <- which(y > 1)
N1 <- length(Id1); N2 <- length(Id2)
Lambda <- c(N1/n, N2/n); MU <- c(mean(y[Id1]), mean(y[Id2]))
Psi <- rep(NA, n); per1 <- c(1,2); per2 <- c(2,1)
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
      if(Psi[i] == 1){Psi[i] <- 2
      }else{Psi[i] <- 1}
    }
  }
  PostPsi[s, ] <- Psi; PostMu[s,] <- MU
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
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
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
  Lambda <- sort(MCMCpack::rdirichlet(1, c(a01 + N1, a02 + N2)), decreasing = TRUE)
  Lambda1 <- Lambda[1]; Lambda2 <- Lambda[2]
  PostLambda[s] <- Lambda1 
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
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

########### Dirichlet process mixture: Simulation ###############
rm(list = ls())
set.seed(010101)
# Simulate data from a 2-component mixture model
N <- 1000; x <- rnorm(N); z <- rbinom(N, 1, 0.5)
y <- ifelse(z == 0, rnorm(N, 2 + 1.5*x, 1), rnorm(N, -1 + 0.5*x, 0.8))
X <- cbind(1, x)
k <- 2
data <- data.frame(y, x); Reg <- lm(y ~ x)
SumReg <- summary(Reg)
# Hyperparameters
a0 <- 0.001; d0 <- 0.001
b0 <- rep(0, k); B0 <- diag(k)
B0i <- solve(B0)
a <- 0.1; b <- 0.1
# MCMC parameters
mcmc <- 5000; burnin <- 1000
tot <- mcmc + burnin; thin <- 2
# Gibbs sampling functions
PostSig2 <- function(Xh, yh){
  Nh <- length(yh)
  yh <- matrix(yh, Nh, 1)
  if(Nh == 1){
    Xh <- matrix(Xh, k, 1)
    Bn <- solve(Xh%*%t(Xh) + B0i)
    bn <- Bn%*%(B0i%*%b0 + Xh%*%yh)
  }else{
    Xh <- matrix(Xh, Nh, k)
    Bn <- solve(t(Xh)%*%Xh + B0i)
    bn <- Bn%*%(B0i%*%b0 + t(Xh)%*%yh)
  }
  Bni <- solve(Bn)
  an <- a0 + Nh
  dn <- d0 + t(yh)%*%yh + t(b0)%*%B0i%*%b0 - t(bn)%*%Bni%*%bn 
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBeta <- function(sig2h, Xh, yh){
  Nh <- length(yh)
  yh <- matrix(yh, Nh, 1)
  if(Nh == 1){
    Xh <- matrix(Xh, k, 1)
    Bn <- solve(Xh%*%t(Xh) + B0i)
    bn <- Bn%*%(B0i%*%b0 + Xh%*%yh)
  }else{
    Xh <- matrix(Xh, Nh, k)
    Bn <- solve(t(Xh)%*%Xh + B0i)
    bn <- Bn%*%(B0i%*%b0 + t(Xh)%*%yh)
  }
  Beta <- MASS::mvrnorm(1, bn, sig2h*Bn)
  return(Beta)
}
PostAlpha <- function(s, alpha){
  H <- length(unique(s))
  psi <- rbeta(1, alpha + 1, N)
  pi.ratio <- (a + H - 1) / (N * (b - log(psi)))
  pi <- pi.ratio / (1 + pi.ratio)
  components <- sample(1:2, prob = c(pi, (1 - pi)), size = 1)
  cs <- c(a + H, a + H - 1)
  ds <- b - log(psi)
  alpha <- rgamma(1, cs[components], ds)
  return(alpha)
}
LogMarLikLM <- function(xh, yh){
  xh <- matrix(xh, k, 1)
  Bn <- solve(xh%*%t(xh) + B0i)
  Bni <- solve(Bn)
  bn <- Bn%*%(B0i%*%b0 + xh%*%yh)
  an <- a0 + 1
  dn <- d0 + yh^2 + t(b0)%*%B0i%*%b0 - t(bn)%*%Bni%*%bn 
  # Log marginal likelihood
  logpy <- (1/2)*log(1/pi)+(a0/2)*log(d0)-(an/2)*log(dn) + 0.5*log(det(Bn)/det(B0)) + lgamma(an/2)-lgamma(a0/2)
  return(logpy)
}
PostS <- function(BETA, SIGMA, Alpha, s, i){
  Nl <- table(s[-i]); H <- length(Nl)
  qh <- sapply(1:H, function(h){(Nl[h]/(N+Alpha-1))*dnorm(y[i], mean = t(X[i,])%*%BETA[,h], sd = SIGMA[h])})
  q0 <- (Alpha/(N+Alpha-1))*exp(LogMarLikLM(xh = X[i,], yh = y[i]))
  # qh <- c(qh, q0)
  qh <- c(q0, qh)
  # si <- sample(1:length(qh), 1, prob = qh)
  Clust <- as.numeric(names(Nl))
  # si <- sample(1:length(qh), 1, prob = qh)
  si <- sample(c(0, Clust), 1, prob = qh)
  # if(si == (H + 1)){
  if(si == 0){
    si <- Clust[H] + 1
    Sig2New <- PostSig2(Xh = X[i,], yh = y[i])
    SIGMA <- c(SIGMA, Sig2New^0.5)
    BetaNew <- PostBeta(sig2h = Sig2New, Xh = X[i,], yh = y[i])
    BETA <- cbind(BETA, BetaNew)
  }
  else {
    si == si
  }
  # return(si)
  return(list(si = si, BETA = BETA, SIGMA = SIGMA))
}
PostBetas <- list(); PostSigma <- list()
Posts <- matrix(0, tot, N); PostAlphas <- rep(0, tot)
S <- sample(1:3, N, replace = T, prob = c(0.5, 0.3, 0.2))
BETA <- cbind(Reg$coefficients, Reg$coefficients, Reg$coefficients)
SIGMA <- rep(SumReg$sigma, 3)
Alpha <- rgamma(1, a, b)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  for(i in 1:N){
    Rests <- PostS(BETA = BETA, SIGMA = SIGMA, Alpha = Alpha, s = S, i = i)
    S[i] <- Rests$si
    BETA <- Rests$BETA; SIGMA <- Rests$SIGMA
  }
  sFreq <- table(S)
  lt <- 1
  for(li in as.numeric(names(sFreq))){
    Index <- which(S == li)
    if(li == lt){
      S[Index] <- li
    } else {
      S[Index] <- lt
    }
    lt <- lt + 1
  }
  Alpha <- PostAlpha(s = S, alpha = Alpha)
  Nl <- table(S); H <- length(Nl)
  SIGMA <- rep(NA, H)
  BETA <- matrix(NA, k, H)
  l <- 1
  for(h in unique(S)){
    Idh <- which(S == h)
    SIGMA[l] <- (PostSig2(Xh = X[Idh, ], yh = y[Idh]))^0.5
    BETA[,l] <- PostBeta(sig2h = SIGMA[l]^2, Xh = X[Idh, ], yh = y[Idh])
    l <- l + 1
  }
  PostBetas[[s]] <- BETA
  PostSigma[[s]] <- SIGMA
  Posts[s, ] <- S
  PostAlphas[s] <- Alpha
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorS<- Posts[keep,]
Clusters <- sapply(1:length(keep), function(i){length(table(PosteriorS[i,]))})
table(Clusters)
NClus <- 3
sapply(1:length(keep), function(i){print(table(PosteriorS[i,]))})
PosteriorSIGMA <- matrix(NA, length(keep), NClus)
l <- 1
for (s in keep){
  PosteriorSIGMA[l,] <- PostSigma[[s]][1:NClus]
  l <- l + 1
}
summary(coda::mcmc(PosteriorSIGMA))
plot(coda::mcmc(PosteriorSIGMA))
PosteriorBeta1 <- matrix(NA, length(keep), k)
j <- 1
for(s in keep){
  PosteriorBeta1[j,] <- PostBetas[[s]][,1]
  j <- j + 1
}
print(summary(coda::mcmc(PosteriorBeta1)))
plot(coda::mcmc(PosteriorBeta1))

PosteriorBeta2 <- matrix(NA, length(keep), k)
j <- 1
for(s in keep){
  PosteriorBeta2[j,] <- PostBetas[[s]][,2]
  j <- j + 1
}
print(summary(coda::mcmc(PosteriorBeta2)))
plot(coda::mcmc(PosteriorBeta2))

PosteriorBeta3 <- matrix(NA, length(keep), k)
j <- 1
for(s in keep){
  PosteriorBeta3[j,] <- PostBetas[[s]][,3]
  j <- j + 1
}
print(summary(coda::mcmc(PosteriorBeta3)))
plot(coda::mcmc(PosteriorBeta3))

########### Dirichlet process mixture: Application (Marijuana consumption in Colombia) ###############
rm(list = ls()); set.seed(010101)
Data <- read.csv("https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/MarijuanaColombia.csv")
attach(Data)
y <- LogMarijuana
X <- as.matrix(cbind(1, Data[,-1]))
Reg <- lm(y ~ X - 1)
SumReg <- summary(Reg)
SumReg
k <- dim(X)[2]
N <- dim(X)[1]
# Plot
library(ggplot2)
ggplot(Data, aes(x = LogMarijuana)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(title = "Density Plot: Marijuana (log) monthly consumption in Colombia", x = "y", y = "Density") +
  theme_minimal()

# Hyperparameters
a0 <- 0.001; d0 <- 0.001
b0 <- rep(0, k); B0 <- diag(k)
B0i <- solve(B0)
a <- 0.1; b <- 0.1
# MCMC parameters
mcmc <- 5000; burnin <- 1000
tot <- mcmc + burnin; thin <- 2
# Gibbs sampling functions
PostSig2 <- function(Xh, yh){
  Nh <- length(yh)
  yh <- matrix(yh, Nh, 1)
  if(Nh == 1){
    Xh <- matrix(Xh, k, 1)
    Bn <- solve(Xh%*%t(Xh) + B0i)
    bn <- Bn%*%(B0i%*%b0 + Xh%*%yh)
  }else{
    Xh <- matrix(Xh, Nh, k)
    Bn <- solve(t(Xh)%*%Xh + B0i)
    bn <- Bn%*%(B0i%*%b0 + t(Xh)%*%yh)
  }
  Bni <- solve(Bn)
  an <- a0 + Nh
  dn <- d0 + t(yh)%*%yh + t(b0)%*%B0i%*%b0 - t(bn)%*%Bni%*%bn 
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBeta <- function(sig2h, Xh, yh){
  Nh <- length(yh)
  yh <- matrix(yh, Nh, 1)
  if(Nh == 1){
    Xh <- matrix(Xh, k, 1)
    Bn <- solve(Xh%*%t(Xh) + B0i)
    bn <- Bn%*%(B0i%*%b0 + Xh%*%yh)
  }else{
    Xh <- matrix(Xh, Nh, k)
    Bn <- solve(t(Xh)%*%Xh + B0i)
    bn <- Bn%*%(B0i%*%b0 + t(Xh)%*%yh)
  }
  Beta <- MASS::mvrnorm(1, bn, sig2h*Bn)
  return(Beta)
}
PostAlpha <- function(s, alpha){
  H <- length(unique(s))
  psi <- rbeta(1, alpha + 1, N)
  pi.ratio <- (a + H - 1) / (N * (b - log(psi)))
  pi <- pi.ratio / (1 + pi.ratio)
  components <- sample(1:2, prob = c(pi, (1 - pi)), size = 1)
  cs <- c(a + H, a + H - 1)
  ds <- b - log(psi)
  alpha <- rgamma(1, cs[components], ds)
  return(alpha)
}
LogMarLikLM <- function(xh, yh){
  xh <- matrix(xh, k, 1)
  Bn <- solve(xh%*%t(xh) + B0i)
  Bni <- solve(Bn)
  bn <- Bn%*%(B0i%*%b0 + xh%*%yh)
  an <- a0 + 1
  dn <- d0 + yh^2 + t(b0)%*%B0i%*%b0 - t(bn)%*%Bni%*%bn 
  # Log marginal likelihood
  logpy <- (1/2)*log(1/pi)+(a0/2)*log(d0)-(an/2)*log(dn) + 0.5*log(det(Bn)/det(B0)) + lgamma(an/2)-lgamma(a0/2)
  return(logpy)
}
PostS <- function(BETA, SIGMA, Alpha, s, i){
  Nl <- table(s[-i]); H <- length(Nl)
  qh <- sapply(1:H, function(h){(Nl[h]/(N+Alpha-1))*dnorm(y[i], mean = t(X[i,])%*%BETA[,h], sd = SIGMA[h])})
  q0 <- (Alpha/(N+Alpha-1))*exp(LogMarLikLM(xh = X[i,], yh = y[i]))
  qh <- c(q0, qh)
  Clust <- as.numeric(names(Nl))
  si <- sample(c(0, Clust), 1, prob = qh)
  if(si == 0){
    si <- Clust[H] + 1
    Sig2New <- PostSig2(Xh = X[i,], yh = y[i])
    SIGMA <- c(SIGMA, Sig2New^0.5)
    BetaNew <- PostBeta(sig2h = Sig2New, Xh = X[i,], yh = y[i])
    BETA <- cbind(BETA, BetaNew)
  }else{si == si
  }
  return(list(si = si, BETA = BETA, SIGMA = SIGMA))
}
PostBetas <- list(); PostSigma <- list()
Posts <- matrix(0, tot, N); PostAlphas <- rep(0, tot)
S <- sample(1:3, N, replace = T, prob = c(0.5, 0.3, 0.2))
BETA <- cbind(Reg$coefficients, Reg$coefficients, Reg$coefficients)
SIGMA <- rep(SumReg$sigma, 3)
Alpha <- rgamma(1, a, b)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  for(i in 1:N){
    Rests <- PostS(BETA = BETA, SIGMA = SIGMA, Alpha = Alpha, s = S, i = i)
    S[i] <- Rests$si
    BETA <- Rests$BETA; SIGMA <- Rests$SIGMA
  }
  sFreq <- table(S)
  lt <- 1
  for(li in as.numeric(names(sFreq))){
    Index <- which(S == li)
    if(li == lt){S[Index] <- li
    } else {S[Index] <- lt
    }
    lt <- lt + 1
  }
  Alpha <- PostAlpha(s = S, alpha = Alpha)
  Nl <- table(S); H <- length(Nl)
  SIGMA <- rep(NA, H)
  BETA <- matrix(NA, k, H)
  l <- 1
  for(h in unique(S)){
    Idh <- which(S == h)
    SIGMA[l] <- (PostSig2(Xh = X[Idh, ], yh = y[Idh]))^0.5
    BETA[,l] <- PostBeta(sig2h = SIGMA[l]^2, Xh = X[Idh, ], yh = y[Idh])
    l <- l + 1
  }
  PostBetas[[s]] <- BETA
  PostSigma[[s]] <- SIGMA
  Posts[s, ] <- S
  PostAlphas[s] <- Alpha
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
ResultsDPMmarijuana <- list(Clusters = Posts, Location = PostBetas, Scales = PostSigma)
save(ResultsDPMmarijuana, file = "ResultsDPMmarijuana.RData")
# load(file = "ResultsDPMmarijuana.RData")
keep <- seq((burnin+1), tot, thin)
PosteriorS<- Posts[keep,]
Clusters <- sapply(1:length(keep), function(i){length(table(PosteriorS[i,]))})
table(Clusters)
NClus <- 3
sapply(1:length(keep), function(i){print(table(PosteriorS[i,]))})
PosteriorSIGMA <- matrix(NA, length(keep), NClus)
l <- 1
for (s in keep){
  PosteriorSIGMA[l,] <- PostSigma[[s]][1:NClus]
  l <- l + 1
}
summary(coda::mcmc(PosteriorSIGMA))
plot(coda::mcmc(PosteriorSIGMA))
PosteriorBeta1 <- matrix(NA, length(keep), k)
j <- 1
for(s in keep){
  PosteriorBeta1[j,] <- PostBetas[[s]][,1]
  j <- j + 1
}
colnames(PosteriorBeta1) <- c("Ct", names(Data[,-1]))
HDI <- HDInterval::hdi(PosteriorBeta1, credMass = 0.95)
HDI
print(summary(coda::mcmc(PosteriorBeta1)))
plot(coda::mcmc(PosteriorBeta1))

PosteriorBeta2 <- matrix(NA, length(keep), k)
j <- 1
for(s in keep){
  PosteriorBeta2[j,] <- PostBetas[[s]][,2]
  j <- j + 1
}
colnames(PosteriorBeta2) <- c("Ct", names(Data[,-1]))
HDI <- HDInterval::hdi(PosteriorBeta2, credMass = 0.95)
HDI
print(summary(coda::mcmc(PosteriorBeta2)))
plot(coda::mcmc(PosteriorBeta2))

PosteriorBeta3 <- matrix(NA, length(keep), k)
j <- 1
for(s in keep){
  PosteriorBeta3[j,] <- PostBetas[[s]][,3]
  j <- j + 1
}
colnames(PosteriorBeta3) <- c("Ct", names(Data[,-1]))
HDI <- HDInterval::hdi(PosteriorBeta3, credMass = 0.95)
HDI
print(summary(coda::mcmc(PosteriorBeta3)))
plot(coda::mcmc(PosteriorBeta3))

DataElast <- data.frame(Marijuana = PosteriorBeta1[,8], Cocaine = PosteriorBeta1[,9],
                        Crack = PosteriorBeta1[,10]) 

dens1 <- ggplot(DataElast, aes(x = Marijuana)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(x = "Marijuana", y = "Density") +
  theme_minimal()

dens2 <- ggplot(DataElast, aes(x = Cocaine)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(x = "Cocaine", y = "Density") +
  theme_minimal()

dens3 <- ggplot(DataElast, aes(x = Crack)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(x = "Crack", y = "Density") +
  theme_minimal()

library(ggpubr)
ggarrange(dens1, dens2, dens3,
          labels = c(
            "A", "B", "C"
          ),
          ncol = 3, nrow = 1,
          legend = "bottom",
          common.legend = TRUE
)
############## Basis functions: Estimates using brms package #################
SplineOwn <- function(x, knots, delta){
  if(knots[1] <= x & x < knots[2]){
    u <- (x - knots[1])/delta
    b <- u^3/6
  }else{
    if(knots[2] <= x & x < knots[3]){
      u <- (x - knots[2])/delta
      b <- (1/6)*(1 + 3*u + 3*u^2 - 3*u^3)
    }else{
      if(knots[3] <= x & x < knots[4]){
        u <- (x - knots[3])/delta
        b <- (1/6)*(4 - 6*u^2 + 3*u^3)
      }else{
        if(knots[4] <= x & x < knots[5]){
          u <- (x - knots[4])/delta
          b <- (1/6)*(1 - 3*u + 3*u^2 - u^3)
        }else{
          b <- 0
        }
      }
    }
  }
  return(b)
}
delta <- 1.5
knotsA <- seq(2, 8, delta)
xA <- seq(2, 8, 0.1)
Ens <- sapply(xA, function(xi) {SplineOwn(xi, knots = knotsA, delta = delta)})
plot(xA, Ens, xlab = "x", ylab = "B-spline", main = "Cubic B-spline comparison: own function vs bs")
require(splines)
BSfunc <- bs(xA, knots = knotsA, degree = 3)
lines(xA, BSfunc[,4], col = "red")

# Cubic B-spline function from scratch
cubic_bspline <- function(x, knots, degree = 3) {
  extended_knots <- c(rep(knots[1], degree), knots, rep(knots[length(knots)], degree))
  num_basis <- length(knots) + degree - 1  
  basis_matrix <- matrix(0, nrow = length(x), ncol = num_basis)
  # Function to compute B-spline basis recursively
  b_spline_basis <- function(x, degree, i, knots) {
    if (degree == 0) {
      return(ifelse(x >= knots[i] & x < knots[i + 1], 1, 0))
    } else {
      left_num <- (x - knots[i])
      left_den <- (knots[i + degree] - knots[i])
      left <- ifelse(left_den != 0, (left_num / left_den) * b_spline_basis(x, degree - 1, i, knots), 0)
      
      right_num <- (knots[i + degree + 1] - x)
      right_den <- (knots[i + degree + 1] - knots[i + 1])
      right <- ifelse(right_den != 0, (right_num / right_den) * b_spline_basis(x, degree - 1, i + 1, knots), 0)
      
      return(left + right)
    }
  }
  
  for (i in 1:num_basis) {
    basis_matrix[, i] <- sapply(x, function(xi) b_spline_basis(xi, degree, i, extended_knots))
  }
  if(x[length(x)] == knots[length(knots)]){
    basis_matrix[length(x), num_basis] <- 1
  }
  return(basis_matrix)
}

delta <- 1.5
knotsA <- seq(2, 8, delta)
xA <- seq(2, 8, 0.1)
basis_matrix <- cubic_bspline(xA, knots = knotsA, degree = 3)
library(splines)
bs_matrix <- bs(xA, knots = knotsA[-c(1, length(knotsA))], degree = 3, intercept = TRUE, Boundary.knots = range(knotsA))
par(mfrow = c(1,2))
matplot(xA, basis_matrix, type = "l", lty = 1, col = rainbow(ncol(basis_matrix)), ylab = "B-spline Basis", xlab = "x", main = "Own function")
matplot(xA, bs_matrix_matrix, type = "l", lty = 1, col = rainbow(ncol(basis_matrix)), ylab = "B-spline Basis", xlab = "x", main = "bs function")

########### Simulation: B-splines ###########
rm(list = ls())
library(ggplot2); library(splines)
# Data generation (same as your original code)
set.seed(010101)
x <- seq(0, 1, 0.001)
ysignal <- 0.4 + 0.25*sin(8*x - 5) + 0.4*exp(-16*(4*x - 2.5)^2)
sig <- 0.15
e <- rnorm(length(ysignal), 0, sd = sig)
y <- ysignal + e
N <- 100
ids <- sort(sample(1:length(ysignal), N))
xobs <- x[ids]
yobs <- y[ids]
knots <- seq(0, 1, 0.25)
BS <- bs(xobs, knots = knots, degree = 3, Boundary.knots = range(x), intercept = FALSE)
# Splines
Spline1 <- 0.56 + BS[,-c(7:8)] %*% rnorm(6, 0, 0.35)
Spline2 <- 0.56 + BS[,-c(7:8)] %*% rnorm(6, 0, 0.35)
Spline3 <- 0.56 + BS[,-c(7:8)] %*% rnorm(6, 0, 0.35)
Spline4 <- 0.56 + BS[,-c(7:8)] %*% rnorm(6, 0, 0.35)
# Create data frames for the true signal, observed data, and Splines
data_true_signal <- data.frame(x = x, y = ysignal, Type = "True Signal")
data_obs <- data.frame(x = xobs, y = yobs, Type = "Observed Data")
# Create separate data frames for each Spline
data_Spline1 <- data.frame(x = xobs, y = Spline1, Type = "Spline 1")
data_Spline2 <- data.frame(x = xobs, y = Spline2, Type = "Spline 2")
data_Spline3 <- data.frame(x = xobs, y = Spline3, Type = "Spline 3")
data_Spline4 <- data.frame(x = xobs, y = Spline4, Type = "Spline 4")

# Combine all data into one data frame for ggplot
data <- rbind(data_true_signal, data_obs, data_Spline1, data_Spline2, data_Spline3, data_Spline4)

# Create the plot using ggplot2
ggplot(data, aes(x = x, y = y)) +
  geom_line(data = subset(data, Type == "True Signal"), aes(color = "True Signal"), linewidth = 1) +
  geom_point(data = subset(data, Type == "Observed Data"), aes(color = "Observed Data"), shape = 16) +
  geom_line(data = subset(data, Type == "Spline 1"), aes(color = "Splines"), linewidth = 1, linetype = "solid") +
  geom_line(data = subset(data, Type == "Spline 2"), aes(color = "Splines"), linewidth = 1, linetype = "solid") +
  geom_line(data = subset(data, Type == "Spline 3"), aes(color = "Splines"), linewidth = 1, linetype = "solid") +
  geom_line(data = subset(data, Type == "Spline 4"), aes(color = "Splines"), linewidth = 1, linetype = "solid") +
  scale_color_manual(values = c("True Signal" = "black", 
                                "Observed Data" = "red", 
                                "Splines" = "blue")) +
  labs(y = "y", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "top")

######## Efects of different knots #########
rm(list = ls())
library(ggplot2)
library(splines)

# Data generation
set.seed(10101)
x <- seq(0, 1, 0.001)
ysignal <- 0.4 + 0.25*sin(8*x - 5) + 0.4*exp(-16*(4*x - 2.5)^2)
sig <- 0.15
e <- rnorm(length(ysignal), 0, sd = sig)
y <- ysignal + e
N <- 100
ids <- sort(sample(1:length(ysignal), N))
xobs <- x[ids]
yobs <- y[ids]

# Generate Fits with different knot placements
knots_list <- list(seq(0, 1, 0.33), seq(0, 1, 0.25), seq(0, 1, 0.2), seq(0, 1, 0.1))
Fits <- list()

for (i in 1:4) {
  BS <- bs(xobs, knots = knots_list[[i]], degree = 3, Boundary.knots = range(x), intercept = FALSE)
  fm <- lm(yobs ~ BS)
  Fits[[i]] <- predict(fm)
}

# Create data frames
data_true_signal <- data.frame(x = x, y = ysignal, Type = "True Signal")
data_obs <- data.frame(x = xobs, y = yobs, Type = "Observed Data")
data_preds <- data.frame(
  x = rep(xobs, 4),
  y = c(Fits[[1]], Fits[[2]], Fits[[3]], Fits[[4]]),
  Type = rep(c("Fit 1", "Fit 2", "Fit 3", "Fit 4"), each = length(xobs))
)

# Combine data
data <- rbind(data_true_signal, data_obs, data_preds)

# Create ggplot
ggplot(data, aes(x = x, y = y, color = Type)) +
  geom_line(data = subset(data, Type == "True Signal"), linewidth = 1) +
  geom_point(data = subset(data, Type == "Observed Data"), shape = 16, size = 2) +
  geom_line(data = subset(data, grepl("Fit", Type)), linewidth = 1, linetype = "solid") +
  scale_color_manual(values = c("True Signal" = "black", 
                                "Observed Data" = "red", 
                                "Fit 1" = "blue", 
                                "Fit 2" = "green", 
                                "Fit 3" = "orange", 
                                "Fit 4" = "purple")) +
  labs(y = "y", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "top")

###### Splines using Stan ########
# library(brms)
# library(ggplot2)
# library(dplyr)
# 
# # Simulate data
# set.seed(123)
# n <- 200
# x <- runif(n, -2, 2)  # Predictor with non-linear effect
# z1 <- rnorm(n)        # Additional regressor 1
# z2 <- rnorm(n)        # Additional regressor 2
# epsilon <- rnorm(n, 0, 0.5)  # Noise term
# 
# # True nonlinear function for x
# f_x <- sin(2 * pi * x) 
# 
# # Generate response variable
# y <- f_x + 0.5 * z1 - 0.3 * z2 + 0.2 * z1 * z2 + epsilon
# 
# # Combine into a dataframe
# mydata <- data.frame(y, x, z1, z2)
# 
# # Fit Bayesian regression using brms
# fit <- brm(y ~ s(x) + z1 + z2 + z1:z2, 
#            data = mydata, 
#            family = gaussian(), 
#            chains = 4, iter = 2000, warmup = 1000, 
#            control = list(adapt_delta = 0.95))
# 
# prior_summary(fit) # Summary of priors
# 
# # Summary of the model
# summary(fit)
# 
# # Extract fitted smooth effect of x
# smooth_pred <- conditional_smooths(fit)
# smooth_data <- as.data.frame(smooth_pred[["mu: s(x)"]])
# 
# # Create the population curve
# pop_curve <- data.frame(x = seq(-2, 2, length.out = 200)) %>%
#   mutate(true_fx = sin(2 * pi * x))
# 
# # Plot estimated smooth effect with population curve
# ggplot() +
#   geom_ribbon(data = smooth_data, aes(x = x, ymin = lower__, ymax = upper__, fill = "Credible Interval"), 
#               alpha = 0.2) +  # Posterior uncertainty
#   geom_line(data = smooth_data, aes(x = x, y = estimate__, color = "Estimated Smooth"), linewidth = 1.2) +  # Estimated smooth
#   geom_line(data = pop_curve, aes(x = x, y = true_fx, color = "True Function"), linetype = "dashed", size = 1.2) +  # True function
#   scale_color_manual(name = "Lines",
#                      values = c("Estimated Smooth" = "blue", "True Function" = "red")) +
#   scale_fill_manual(name = "Shaded Area", values = c("Credible Interval" = "blue")) +
#   labs(title = "Estimated Smooth Effect of x vs. True Function",
#        x = "x", y = "f(x)",
#        caption = "Blue: Estimated | Red Dashed: True Function") +
#   theme_minimal()

