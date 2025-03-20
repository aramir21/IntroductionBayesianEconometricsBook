########################## Simulation exercise: Semi-parametric model ########################## 
rm(list = ls())
set.seed(010101)
# Simulate data from a 2-component mixture model
n <- 1000
x1 <- rnorm(n); x2 <- rnorm(n)
B <- c(0.5, -1.2)
X <- cbind(x1, x2)
z <- rbinom(n, 1, 0.3)  # Latent class indicator
mu <- ifelse(z == 0, rnorm(n, -0.5, 0.5), rnorm(n, 1, 0.8))
plot(density(mu))
y <- X%*%B + mu
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
Reg <- lm(y ~ X)
Beta <- Reg$coefficients[2:3]
Res <- Reg$residuals
plot(density(Res))
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

PosteriorSIGMA <- matrix(NA, length(keep), 2)
PosteriorMU <- matrix(NA, length(keep), 2)
PosteriorLAMBDA <- matrix(NA, length(keep), 2)
l <- 1
for (s in keep){
  PosteriorSIGMA[l,] <- PostSigma2[[s]][1:2]
  PosteriorMU[l,] <- PostMu[[s]][1:2]
  PosteriorLAMBDA[l,] <- PostLambda[[s]][1:2]
  l <- l + 1
}

summary(coda::mcmc(PosteriorSIGMA))
summary(coda::mcmc(PosteriorMU))
summary(coda::mcmc(PosteriorLAMBDA))

# Posterior density estimate #
mus <- seq(-3, 4, 0.01)
DensityEst <- function(par, eval){
  lambda1 <- par[1]; lambda2 <- par[2]
  mu1 <- par[3]; mu2 <- par[4]
  sd1 <- par[5]; sd2 <- par[6]
  Dens <- lambda1 * dnorm(eval, mu1, sd1) + lambda2 * dnorm(eval, mu2, sd2)
  return(Dens)
}
par <- cbind(PosteriorLAMBDA, PosteriorMU, PosteriorSIGMA^0.5)
DensEval <- matrix(NA, length(keep), length(mus))
DensPop <- rep(NA, length(mus))
for(r in 1:length(mus)){
  for(l in 1:length(keep)){
    DensEval[l, r] <- DensityEst(par = par[l,], eval = mus[r])
  }
  DensPop[r] <- DensityEst(par = c(0.7, 0.3, -0.5, 1, 0.5, 0.8), eval = mus[r]) 
}
library(dplyr)
library(ggplot2)
require(latex2exp)
DataDens <- tibble(t = mus,
                   Pop = DensPop,
                   lower = apply(DensEval, 2, quantile, probs = 0.025),
                   upper = apply(DensEval, 2, quantile, probs = 0.975),
                   meanT = colMeans(DensEval))
plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1, fill = "lightblue") +
    geom_line(aes(y = Pop, color = "True Value"), linewidth = 0.5) +
    geom_line(aes(y = meanT, color = "Estimate"), linewidth = 0.5) +
    scale_color_manual(values = c("True Value" = "black", "Estimate" = "blue")) +
    xlab(TeX("$//mu$")) + 
    ylab("Density") +
    labs(title = "Density: Stochastic Errors", color = "") + # Label for legend
    theme_minimal()
  print(p)
}
plot_filtering_estimates(DataDens)

########################## Mixture model: Demand of marijuana ########################## 
rm(list = ls()); set.seed(010101)
Data <- read.csv("https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/MarijuanaColombia.csv")
attach(Data)
y <- LogMarijuana
X <- as.matrix(cbind(1, Data[,-1]))
Reg <- lm(y ~ X - 1)
summary(Reg)
k <- dim(X)[2]
N <- dim(X)[1]
# Plot
library(ggplot2)
ggplot(Data, aes(x = LogMarijuana)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(title = "Density Plot", x = "y", y = "Density") +
  theme_minimal()

# Hyperparameters
d0 <- 0.001
a0 <- 0.001
b0 <- rep(0, k)
B0 <- diag(k)
B0i <- solve(B0)
H <- 5
a0h <- rep(1/H, H)
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
  if(length(yh) == 1){
    Xh <- matrix(Xh, k, 1)
    Bn <- solve(B0i + sig2h^(-1)*Xh%*%t(Xh))
    bn <- Bn%*%(B0i%*%b0 + sig2h^(-1)*Xh%*%yh)
  }else{
    Bn <- solve(B0i + sig2h^(-1)*t(Xh)%*%Xh)
    bn <- Bn%*%(B0i%*%b0 + sig2h^(-1)*t(Xh)%*%yh)
  }
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}

PostBetas <- list()
PostSigma2 <- list()
PostLambda <- list()
PostPsi <- matrix(0, tot, N)
Psi <- sample(1:H, N, replace = TRUE, prob = a0h)
Reg <- lm(y ~ X -1)
Beta <- matrix(rep(Reg$coefficients, H), k, H) 
Hs <- H
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Sig2 <- rep(NA, Hs)
  for(h in 1:Hs){
    idh <- which(Psi == h)
    Sig2[h] <- PostSig2(Betah = Beta[, h], Xh = X[idh,], yh = y[idh])
  }
  Hs <- length(Sig2)
  Beta <- matrix(NA, k, Hs)
  for(h in 1:Hs){
    idh <- which(Psi == h)
    Beta[, h] <- PostBeta(sig2h = Sig2[h], Xh = X[idh,], yh = y[idh])
  }
  Lambda <- sort(MCMCpack::rdirichlet(1, a0h[1:Hs] + table(Psi)), decreasing = TRUE)
  for(i in 1:N){
    lambdai <- NULL
    for(h in 1:Hs){
      idh <- which(Psi == h)
      lambdaih <- Lambda[h]*dnorm(y[i], X[i,]%*%Beta[,h], Sig2[h]^0.5)
      lambdai <- c(lambdai, lambdaih)
    }
    Psi[i] <- sample(1:Hs, 1, prob = lambdai)
  }
  Hs <- length(table(Psi))
  PostBetas[[s]] <- Beta
  PostSigma2[[s]] <- Sig2
  PostLambda[[s]] <- Lambda
  PostPsi[s, ] <- Psi
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorPsi <- PostPsi[keep,]
Clusters <- sapply(1:length(keep), function(i){length(table(PosteriorPsi[i,]))})
NClus <- 4
# Clusters <- sapply(1:length(keep), function(i){print(table(PosteriorPsi[i,]))})
PosteriorSIGMA <- matrix(NA, length(keep), NClus)
PosteriorBETA <- array(NA, c(k,NClus,length(keep)))
PosteriorLAMBDA <- matrix(NA, length(keep), NClus)
l <- 1
for (s in keep){
  PosteriorSIGMA[l,] <- PostSigma2[[s]][1:NClus]
  PosteriorLAMBDA[l,] <- PostLambda[[s]][1:NClus]
  PosteriorBETA[,,l] <- PostBetas[[s]]
  l <- l + 1
}

summary(coda::mcmc(PosteriorSIGMA))
plot(coda::mcmc(PosteriorSIGMA))
summary(coda::mcmc(PosteriorLAMBDA))
plot(coda::mcmc(PosteriorLAMBDA))
for(l in 1:NClus){
  PosteriorBeta <- t(PosteriorBETA[,l,])
  colnames(PosteriorBeta) <- c("Ct", names(Data[,-1]))
  print(summary(coda::mcmc(PosteriorBeta)))
}
PosteriorBeta1 <- t(PosteriorBETA[,1,])
colnames(PosteriorBeta1) <- c("Ct", names(Data[,-1]))
plot(coda::mcmc(PosteriorBeta1))

PosteriorBeta2 <- t(PosteriorBETA[,2,])
colnames(PosteriorBeta2) <- c("Ct", names(Data[,-1]))
plot(coda::mcmc(PosteriorBeta2))

########### Dirichlet process mixture: Semi-parametric simulation ###############
rm(list = ls())
set.seed(010101)
# Simulate data from a 2-component mixture model
n <- 1000
x1 <- rnorm(n); x2 <- rnorm(n)
B <- c(0.5, -1.2)
X <- cbind(x1, x2)
z <- rbinom(n, 1, 0.3)  # Latent class indicator
mu <- ifelse(z == 0, rnorm(n, -0.5, 0.5), rnorm(n, 1, 0.8))
y <- X%*%B + mu
k <- dim(X)[2]
N <- dim(X)[1]
# Hyperparameters
a0 <- 0.001; d0 <- 0.001
b0 <- rep(0, k); B0 <- diag(k)
B0i <- solve(B0)
a <- 0.1; b <- 0.1
mu0 <- 0
sig2mu0 <- 10
beta0 <- 0.1
# MCMC parameters
mcmc <- 5000; burnin <- 1000
tot <- mcmc + burnin; thin <- 2
# Gibbs sampling functions
PostSig2 <- function(Xh, yh, Beta){
  Nh <- length(yh)
  if(Nh == 1){
    Xh <- matrix(Xh, k, 1)
    XhB <- t(Xh)%*%Beta
  }else{
    Xh <- matrix(Xh, Nh, k)
    XhB <- Xh%*%Beta
  }
  eh <- matrix(yh, Nh, 1) - XhB
  an <- a0 + Nh
  dn <- d0 + sum((eh - mean(eh))^2) + beta0*Nh/(beta0 + Nh)*(mean(eh) - mu0)^2
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
Postmu <- function(sig2h, Beta, Xh, yh){
  Nh <- length(yh)
  if(Nh == 1){
    Xh <- matrix(Xh, k, 1)
    XhB <- t(Xh)%*%Beta
  }else{
    Xh <- matrix(Xh, Nh, k)
    XhB <- Xh%*%Beta
  }
  eh <- matrix(yh, Nh, 1) - XhB
  sig2mu <- sig2h/(beta0 + Nh)
  mun <- (beta0*mu0 + Nh*mean(eh))/(beta0 + Nh)
  mu <- rnorm(1, mun, sig2mu^0.5)
  return(mu)
}
PostBeta <- function(sig2, mu, X, y, Clust){
  XtX <- matrix(0, k, k)
  Xty <- matrix(0, k, 1)
  Hs <- length(mu)
  for(h in 1:Hs){
    idh <- which(Clust == h)
    if(length(idh) == 1){
      Xh <- matrix(X[idh,], k, 1)
      XtXh <- sig2[h]^(-1)*Xh%*%t(Xh)
      yh <- y[idh]
      Xtyh <- sig2[h]^(-1)*Xh%*%(yh - mu[h])
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
PostAlpha <- function(Clust, alpha){
  H <- length(unique(Clust))
  psi <- rbeta(1, alpha + 1, N)
  pi.ratio <- (a + H - 1) / (N * (b - log(psi)))
  pi <- pi.ratio / (1 + pi.ratio)
  components <- sample(1:2, prob = c(pi, (1 - pi)), size = 1)
  cs <- c(a + H, a + H - 1)
  ds <- b - log(psi)
  alpha <- rgamma(1, cs[components], ds)
  return(alpha)
}
LogMarLikLM <- function(xh, yh, Beta){
  xh <- matrix(xh, k, 1)
  ei <- c(yh - t(xh)%*%Beta)
  betan <- beta0 + 1
  an <- a0 + 1
  dn <- d0 +beta0/(beta0 + 1)*(ei - mu0)^2  
  # Log marginal likelihood
  logpy <- (1/2)*log(1/pi)+(a0/2)*log(d0)-(an/2)*log(dn) + 0.5*log(beta0/betan) + lgamma(an/2)-lgamma(a0/2)
  return(logpy)
}
PostS <- function(Beta, SIGMA, MU, Alpha, Clust, i){
  Nl <- table(Clust[-i]); H <- length(Nl)
  ei <- c(y[i] - t(X[i,])%*%Beta)
  qh <- sapply(1:H, function(h){(Nl[h]/(N+Alpha-1))*dnorm(ei, mean = MU[h], sd = SIGMA[h])})
  q0 <- (Alpha/(N+Alpha-1))*exp(LogMarLikLM(xh = X[i,], yh = y[i], Beta = Beta))
  qh <- c(q0, qh)
  Clust <- as.numeric(names(Nl))
  si <- sample(c(0, Clust), 1, prob = qh)
  if(si == 0){
    si <- Clust[H] + 1
    Sig2New <- PostSig2(Xh = X[i,], yh = y[i], Beta = Beta)
    SIGMA <- c(SIGMA, Sig2New^0.5)
    MuNew <- Postmu(sig2h = Sig2New, Beta = Beta, Xh = X[i,], yh = y[i])
    MU <- cbind(MU, MuNew)
  }else{si == si}
  return(list(si = si, MU = MU, SIGMA = SIGMA))
}
PostBetas <- matrix(NA, tot, k)
PostSigma <- list(); PostMu <- list()
Posts <- matrix(0, tot, N); PostAlphas <- rep(0, tot)
Clust <- sample(1:3, N, replace = T, prob = c(0.5, 0.3, 0.2))
Reg <- lm(y ~ X)
SumReg <- summary(Reg)
BETA <- Reg$coefficients[-1]
SIGMA <- rep(SumReg$sigma, 3)
MU <- rep(Reg$coefficients[1], 3)
Alpha <- rgamma(1, a, b)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  for(i in 1:N){
    Rests <- PostS(Beta = BETA, SIGMA = SIGMA, MU = MU, Alpha = Alpha, Clust = Clust, i = i)
    Clust[i] <- Rests$si
    MU <- Rests$MU; SIGMA <- Rests$SIGMA
  }
  sFreq <- table(Clust)
  lt <- 1
  for(li in as.numeric(names(sFreq))){
    Index <- which(Clust == li)
    if(li == lt){Clust[Index] <- li
    } else {Clust[Index] <- lt
    }
    lt <- lt + 1
  }
  Alpha <- PostAlpha(Clust = Clust, alpha = Alpha)
  Nl <- table(Clust); H <- length(Nl)
  SIGMA <- rep(NA, H)
  MU <- rep(NA, H)
  l <- 1
  for(h in unique(Clust)){
    Idh <- which(Clust == h)
    SIGMA[l] <- (PostSig2(Beta = BETA, Xh = X[Idh, ], yh = y[Idh]))^0.5
    MU[l] <- Postmu(Beta = BETA, sig2h = SIGMA[l]^2, Xh = X[Idh, ], yh = y[Idh])
    l <- l + 1
  }
  BETA <- PostBeta(sig2 = SIGMA, mu = MU, X = X, y = y, Clust = Clust) 
  PostBetas[s,] <- BETA
  PostSigma[[s]] <- SIGMA
  PostMu[[s]] <- MU
  Posts[s, ] <- Clust
  PostAlphas[s] <- Alpha
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq(burnin, tot, thin)
PosteriorBetas <- coda::mcmc(PostBetas[keep,])
summary(PosteriorBetas)
plot(PosteriorBetas)
PosteriorPsi <- Posts[keep,]
Clusters <- sapply(1:length(keep), function(i){length(table(PosteriorPsi[i,]))})
table(Clusters)
PosteriorSIGMA <- matrix(NA, length(keep), 2)
PosteriorMU <- matrix(NA, length(keep), 2)
l <- 1
for (s in keep){
  PosteriorSIGMA[l,] <- PostSigma[[s]][1:2]
  PosteriorMU[l,] <- PostMu[[s]][1:2]
  l <- l + 1
}
summary(coda::mcmc(PosteriorSIGMA))
summary(coda::mcmc(PosteriorMU))
plot(coda::mcmc(PosteriorMU))
DataMU2 <- data.frame(Posterior = PosteriorMU[,2]) 
library(ggplot2)
densMU <- ggplot(DataMU2, aes(x = Posterior)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(x = "Intercept", y = "Density") +
  theme_minimal()
densMU 

############### Dirichlet process solution switching label issue using finite mixtures #######################
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
H <- 3
a0h <- rep(1/H, H)
# MCMC parameters
mcmc <- 10000
burnin <- 2000
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
PostBetas3 <- matrix(0, tot, 2)
PostSigma2 <- matrix(0, tot, H)
PostPsi <- matrix(0, tot, n)
PostLambda <- matrix(0, tot, H)
Reg <- lm(y ~ x)
SumReg <- summary(Reg)
BETA <- cbind(Reg$coefficients, Reg$coefficients, Reg$coefficients)
SIG2 <- rep(SumReg$sigma, H)
X <- cbind(1, x)
Psi <- rep(NA, n)
Lambda1 <- 1/H; Lambda2 <- 1/H; Lambda3 <- 1/H
perutatms <- gtools::permutations(3, 3, 1:H)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  for(i in 1:n){
    lambdai1 <- Lambda1*dnorm(y[i], X[i,]%*%BETA[,1], SIG2[1]^0.5)
    lambdai2 <- Lambda2*dnorm(y[i], X[i,]%*%BETA[,2], SIG2[2]^0.5)
    lambdai3 <- Lambda3*dnorm(y[i], X[i,]%*%BETA[,3], SIG2[3]^0.5)
    Psi[i] <- sample(1:H, 1, prob = c(lambdai1, lambdai2, lambdai3))
  }
  NH <- NULL
  for(h in 1:H){
    Idh <- which(Psi == h)
    Nh <- length(Idh)
    NH <- c(NH, Nh)
    SIG2[h] <- PostSig2(Betah = BETA[,h], Xh = X[Idh, ], yh = y[Idh])
    BETA[,h] <- PostBeta(sig2h = SIG2[h], Xh = X[Idh, ], yh = y[Idh])
  }
  Lambda <- sort(MCMCpack::rdirichlet(1, a0h + NH), decreasing = TRUE)
  IndPer <- sample(1:dim(perutatms)[1], 1)
  Clust <- perutatms[IndPer, ]
  SIG2 <- SIG2[Clust]
  BETA <- BETA[,Clust]
  Lambda <- Lambda[Clust]
  for(i in 1:n){
    if(Psi[i] == 1){
      Psi[i] <- Clust[1]
    }else{
      if(Psi[i] == 2){
        Psi[i] <- Clust[2]
      }else{
        Psi[i] <- Clust[3]
      }
    }
  }
  Lambda1 <- Lambda[1]; Lambda2 <- Lambda[2]; Lambda3 <- Lambda[3]
  PostPsi[s, ] <- Psi
  PostSigma2[s, ] <- SIG2
  PostBetas1[s,] <- BETA[,1]
  PostBetas2[s,] <- BETA[,2]
  PostBetas3[s,] <- BETA[,3]
  PostLambda[s,] <- Lambda 
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
PosteriorBetas3 <- coda::mcmc(PostBetas3[keep,])
summary(PosteriorBetas3)
plot(PosteriorBetas3)
PosteriorSigma2 <- coda::mcmc(PostSigma2[keep,])
summary(PosteriorSigma2)
plot(PosteriorSigma2)
PosteriorLambda <- coda::mcmc(PostLambda[keep,])
summary(PosteriorLambda)
plot(PosteriorLambda)

DataBetas <- data.frame(Cluster1 = as.numeric(PosteriorBetas1[,2]), Cluster2 = as.numeric(PosteriorBetas2[,2]),
                        Cluster3 = as.numeric(PosteriorBetas3[,2])) 
library(ggplot2)
dens1 <- ggplot(DataBetas, aes(x = Cluster1)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(x = "Cluster 1", y = "Density") +
  theme_minimal()

dens2 <- ggplot(DataBetas, aes(x = Cluster2)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(x = "Cluster 2", y = "Density") +
  theme_minimal()

dens3 <- ggplot(DataBetas, aes(x = Cluster3)) +
  geom_density(fill = "blue", alpha = 0.3) +  # Density plot with fill color
  labs(x = "Cluster 3", y = "Density") +
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

########### Splines: Application (Marijuana consumption in Colombia) ###############
rm(list = ls()); 
library(splines); library(ggplot2)
set.seed(010101)
Data <- read.csv("https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/MarijuanaColombia.csv")
attach(Data)
IdOrd <- order(Age) 
y <- LogMarijuana[IdOrd]
Z <- as.matrix(cbind(Data[IdOrd,-c(1, 6, 7)]))
x <- Age[IdOrd] 
knots <- quantile(x, seq(0, 1, 0.05))
BS <- bs(x, knots = knots, degree = 3, Boundary.knots = range(x), intercept = FALSE)
matplot(x, BS, type = "l", lty = 1, col = rainbow(ncol(BS)))
# Function to check if a column is constant
is_constant <- function(col) {
  return(length(unique(col)) == 1)
}
constant_columns <- apply(BS, 2, is_constant)
constant_columns
X <- cbind(BS[,-c(which(constant_columns == TRUE))], Z)
BMAglm <- BMA::bicreg(X, y, strict = FALSE, OR = 50) 
summary(BMAglm)
idReg <- c(2, 23, 24, 27:29) # Relevant regressors, PIP > 0.5
Xnew <- cbind(1, X[, idReg])
N <- dim(Xnew)[1]
k <- dim(Xnew)[2]
# Hyperparameters
d0 <- 0.001
a0 <- 0.001
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
# MCMC parameters
mcmc <- 5000
burnin <- 5000
tot <- mcmc + burnin
thin <- 1
# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
posterior  <- MCMCpack::MCMCregress(y~Xnew-1, b0=b0, B0 = B0i, c0 = a0, d0 = d0, burnin = burnin, mcmc = mcmc, thin = thin)
summary(coda::mcmc(posterior))
# Predict values with 95% credible intervals
xfit <- seq(min(x), max(x), 0.2)
H <- length(xfit)
i <- 675
idfit <- sample(1:N, H)
BSfit <- bs(xfit, knots = knots, degree = 3, Boundary.knots = range(x), intercept = FALSE)
Xfit <- cbind(1, BSfit[,c(2,23)], Z[rep(i, H),c(1,4:6)]) # Relevant regressors, PIP > 0.5
Fit <- matrix(NA, mcmc, H)
# posterior[posterior > 0] <- 0
for(s in 1:mcmc){
  Fit[s,] <- Xfit%*%posterior[s,1:7]
}
# Create a data frame for ggplot
plot_data <- data.frame(x = xfit, fit = colMeans(Fit), liminf = apply(Fit, 2, quantile, 0.025), 
                        limsup = apply(Fit, 2, quantile, 0.975))

ggplot() +
  geom_line(data = plot_data, aes(x, fit), color = "blue", linewidth = 1) +  # Regression line
  geom_ribbon(data = plot_data, aes(x, ymin = liminf, ymax = limsup), fill = "blue", alpha = 0.2) +  # Confidence interval
  labs(title = "B-Spline Regression with 95% Confidence Interval",
       x = "Age",
       y = "Log Marijuana") +
  theme_minimal()

