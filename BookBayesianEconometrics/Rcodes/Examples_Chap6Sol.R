rm(list = ls())
set.seed(010101)
########################## Linear regression: Value of soccer players (heteroskedasticity) ##########################
Data <- read.csv("DataApplications/1ValueFootballPlayers.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(Data)
y <- log(Value) 
# Value: Market value in Euros (2017) of soccer players
# Regressors quantity including intercept
X <- cbind(1, Perf, Age, Age2, NatTeam, Goals, Exp, Exp2)
# Perf: Performance. Perf2: Performance squared. Age: Age; Age: Age squared. 
# NatTeam: Indicator of national team. Goals: Scored goals. Goals2: Scored goals squared
# Exp: Years of experience. Exp2: Years of experience squared. Assists: Number of assists
k <- dim(X)[2]
N <- dim(X)[1]
# Hyperparameters
d0 <- 0.001
a0 <- 0.001
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
v <- 5
# MCMC parameters
mcmc <- 5000
burnin <- 5000
tot <- mcmc + burnin
thin <- 1
# Posterior distributions programming the Gibbs sampling
# Auxiliary parameters
an <- a0 + N
v1n <- v + 1
# Gibbs sampling functions
PostSig2 <- function(Beta, tau){
  dn <- d0 + t(y - X%*%Beta)%*%diag(tau)%*%(y - X%*%Beta)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBeta <- function(sig2, tau){
  Bn <- solve(B0i + sig2^(-1)*t(X)%*%diag(tau)%*%X)
  bn <- Bn%*%(B0i%*%b0 + sig2^(-1)*t(X)%*%diag(tau)%*%y)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostTau <- function(sig2, Beta, i){
  v2n <- v + sig2^(-1)*(y[i]-X[i,]%*%Beta)^2
  taui <- rgamma(1, v1n/2, v2n/2)
  return(taui)
}
PostBetas <- matrix(0, mcmc+burnin, k)
PostSigma2 <- rep(0, mcmc+burnin)
Beta <- rep(0, k)
tau <- rep(1, N)
# create progress bar in case that you want to see iterations progress
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  sig2 <- PostSig2(Beta = Beta, tau = tau)
  PostSigma2[s] <- sig2
  Beta <- PostBeta(sig2 = sig2, tau = tau)
  PostBetas[s,] <- Beta
  tau <- sapply(1:N, function(i){PostTau(sig2 = sig2, Beta = Beta, i)})
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0), "% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
colnames(PosteriorBetas) <- c("Intercept", "Perf", "Age", "Age2", "NatTeam", "Goals", "Exp", "Exp2")
summary(coda::mcmc(PosteriorBetas))
PosteriorSigma2 <- PostSigma2[keep]
summary(coda::mcmc(PosteriorSigma2))
summary(coda::mcmc(exp(PosteriorBetas[,5])-1))

########################## Probit model: Determinants of hospitalization ##########################
rm(list = ls())
set.seed(010101)
Data <- read.csv("DataApplications/2HealthMed.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(Data)
str(Data)
y <- Hosp # Dependent variables
X <- cbind(1, SHI, Female, Age, Age2, Est2, Est3, Fair, Good, Excellent) # Regressors
K <- dim(X)[2] # Number of regressors
N <- dim(X)[1]
# Hyperparameters
b0 <- rep(0, K) # Prio mean
B0 <- diag(K) # Prior precision (inverse of covariance)
B0i <- solve(B0)
mcmc <- 1000; burnin <- 500; thin <- 2; tot <- mcmc + burnin; keep <- seq(burnin, tot, thin)
# Posterior distributions programming the Gibbs sampling
# Auxiliary parameters
XtX <- t(X)%*%X
# Gibbs sampling functions
PostBeta <- function(Yl){
  Bn <- solve(B0i + XtX)
  bn <- Bn%*%(B0i%*%b0 + t(X)%*%Yl)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostYl <- function(Beta, i){
  Ylmean <- X[i,]%*%Beta
  if(y[i] == 1){
    Yli <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = Ylmean, sd = 1)
  }else{
    Yli <- truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = Ylmean, sd = 1)
  }
  return(Yli)
}
PostBetas <- matrix(0, mcmc+burnin, K)
Beta <- rep(0, K)
# create progress bar in case that you want to see iterations progress
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Yl <- sapply(1:N, function(i){PostYl(Beta = Beta, i)})
  Beta <- PostBeta(Yl = Yl)
  PostBetas[s,] <- Beta
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0), "% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
colnames(PosteriorBetas) <- c("Intercept", "SHI", "Female", "Age", "Age2", "Est2", "Est3", "Fair", "Good", "Excellent")
summary(coda::mcmc(PosteriorBetas))

########################## Multinomial probit: Fishing choice ##########################
remove(list = ls())
set.seed(12345)
Data <- read.csv("DataApplications/3Fishing.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(Data)
str(Data)
p <- 4; na <- 2; nd <- 1; N <- dim(Data)[1]
Price <- c(t(Data[,2:5]))
Catch <- c(t(Data[,6:9]))
Xa <- Data[,2:9]
Xd <- matrix(income, N, 1)
X <- bayesm::createX(p = p, na = na, nd = nd, Xa = Xa, Xd = Xd, INT = TRUE, base = p, DIFF = TRUE)
df <- list(y = mode, X = X, p = 4)
# Hyperparameters
k <- dim(X)[2]
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
a0 <- p - 1 + 3
Psi0 <- a0*diag(p-1)
Prior <- list(betabar = b0, A = B0i, nu = a0, V = Psi0)
# MCMC parameters
mcmc <- 100000
thin <- 5
Mcmc <- list(R = mcmc, keep = thin)
Results <- bayesm::rmnpGibbs(Data = df, Prior = Prior, Mcmc = Mcmc)

betatilde <- Results$betadraw / sqrt(Results$sigmadraw[,1])
attributes(betatilde)$class <- "bayesm.mat"
summary(coda::mcmc(betatilde))

sigmadraw <-  Results$sigmadraw / Results$sigmadraw[,1]
attributes(sigmadraw)$class = "bayesm.var"
summary(coda::mcmc(sigmadraw))

########################## Multinomial Logit: Simulation ##########################
remove(list = ls())
set.seed(12345)
# Simulation of data
N<-1000  # Sample Size
B<-c(0.5,0.8,-3)
B1<-c(-2.5,-3.5,0)
B2<-c(1,1,0)
# Alternative specific attributes of choice 1, for instance, price, quality and duration of choice 1
X1<-matrix(cbind(rnorm(N,0,1),rnorm(N,0,1),rnorm(N,0,1)),N,length(B)) 
# Alternative specific attributes of choice 2, for instance, price, quality and duration of choice 2
X2<-matrix(cbind(rnorm(N,0,1),rnorm(N,0,1),rnorm(N,0,1)),N,length(B))
# Alternative specific attributes of choice 3, for instance, price, quality and duration of choice 3
X3<-matrix(cbind(rnorm(N,0,1),rnorm(N,0,1),rnorm(N,0,1)),N,length(B))
X4<-matrix(rnorm(N,1,1),N,1)
V1<-B2[1]+X1%*%B+B1[1]*X4
V2<-B2[2]+X2%*%B+B1[2]*X4
V3<-B2[3]+X3%*%B+B1[3]*X4
suma<-exp(V1)+exp(V2)+exp(V3)
p1<-exp(V1)/suma
p2<-exp(V2)/suma
p3<-exp(V3)/suma
p<-cbind(p1,p2,p3)
y<- apply(p,1, function(x)sample(1:3, 1, prob = x, replace = TRUE))
table(y)
L <- length(table(y))
dat <-data.frame(y,X1[,1],X2[,1],X3[,1],X1[,2],X2[,2],X3[,2],X1[,3],X2[,3],X3[,3],X4)
colnames(dat) <- c("mode","V1.1","V1.2","V1.3","V2.1","V2.2","V2.3","V3.1","V3.2","V3.3","V4")
attach(dat)
LongData <- mlogit::mlogit.data(dat, shape = "wide", varying=2:10, choice = "mode")
Xa <- cbind(LongData$V1, LongData$V2, LongData$V3)
Xa <- cbind(X1[,1],X2[,1],X3[,1],X1[,2],X2[,2],X3[,2],X1[,3],X2[,3],X3[,3])
na <- 3
Xd <- X4
X <- bayesm::createX(p = L, na = na, nd = 1, Xa = Xa, Xd = Xd, base = L)
DataMlogit <- list(y=y, X = X, p = L)
# MCMC parameters
mcmc <- 11000+1
thin <- 5
df <- 6
mcmcpar <- list(R = mcmc, keep = 5, nu = df)
PostBeta <- bayesm::rmnlIndepMetrop(Data = DataMlogit, Mcmc = mcmcpar)
summary(PostBeta[["betadraw"]])

########################## Ordered probit: Simulation ##########################
rm(list = ls())
set.seed(010101)
N <- 1000
x1 <- rnorm(N, 6, 5); x2 <- rgamma(N, shape = 1, scale = 1)
X <- cbind(1, x1, x2)
beta <- c(0.5, -0.25, 0.5); cutoffs <- c(0, 1, 2.5)
e <- rnorm(N,0,1)
y_latent <- X%*%beta + e
y <- rep(0,N)
for (i in 1:N) {
  if (y_latent[i] < cutoffs[1]) {
    y[i] <- 0
  } else {
    if (y_latent[i] >= cutoffs[1] & y_latent[i] < cutoffs[2]) {
      y[i] <- 1
    } else{
      if (y_latent[i] >= cutoffs[2] & y_latent[i] < cutoffs[3]) {
        y[i] <- 2
      } else{
        y[i] <- 3
      }
    }
  }
}
table(y)
# Likelihood function
LogLikOP <- function(param){
  beta_g <- param[1:ncol(X)]
  delta <- param[(ncol(X)+1):(ncol(X) + dplyr::n_distinct(y) - 1)]
  Xbeta <- X%*%beta_g
  logLik <- 0
  for (i in 1:length(y)) {
    if (y[i]==0) {
      logLiki <- log(pnorm(-Xbeta[i]))
    } else if (y[i]==1) {
      logLiki <- log(pnorm(exp(delta[1]) - Xbeta[i]) - pnorm(-Xbeta[i]))
    } else if (y[i]==2){
      logLiki <- log(pnorm(exp(delta[2]) + exp(delta[1]) - Xbeta[i]) - pnorm(exp(delta[1]) - Xbeta[i]))
    } else {
      logLiki <- log(1 - pnorm(exp(delta[2]) + exp(delta[1]) - Xbeta[i]))
    }
    logLik <- logLik + logLiki
  }
  return(-logLik)
}
# ML Estimation
param0 <- rep(0, ncol(X) + n_distinct(y)-2)
mle <- optim(param0, LogLikOP, hessian = T, method = "BFGS")
mle$par
exp(mle$par[length(beta)+1])
exp(mle$par[length(beta)+1])+exp(mle$par[length(beta)+2])
CovarML <- solve(mle$hessian)
# M-H within Gibbs
mhop <- function(param0, G){
  betasamples <- matrix(c(0), nrow = G, ncol = ncol(X))
  betasamples[1,] <- param0[1:ncol(X)]
  tau <- matrix(c(0), nrow = G, ncol = dplyr::n_distinct(y) - 2)
  tau[1,] <- param0[(ncol(X)+1):(ncol(X) + dplyr::n_distinct(y) - 2)]
  yl <- rep(0,length(y))
  ar <- rep(0,G)
  B1 <- solve(t(X)%*%X+solve(B0))
  pb <- winProgressBar(title = "progress bar", min = 0, max = G, width = 300)
  for(g in 2:G){
    bg <- betasamples[g-1,]
    tg <- tau[g-1,]
    #Random walk M-H for delta
    delta_prime <- tg + mvtnorm::rmvnorm(1, mean = rep(0,2), sigma = VarProp)
    alpha <- min(1,(mvtnorm::dmvnorm(delta_prime, mean = d0, sigma = D0)*exp(-LogLikOP(c(bg, delta_prime)) + LogLikOP(c(bg, tg))))/mvtnorm::dmvnorm(tg, mean = d0, sigma = D0))
    if(is.nan(alpha) | is.na(alpha)) {
      alpha <- 0
    }
    #Acceptance step
    u <- runif(1, min = 0, max = 1)
    if(u<=alpha){
      tau[g,] <- delta_prime
      ar[g] <- 1
    }else{
      tau[g,] <- tg
    }
    #Generation of latent variables
    for (i in 1:length(y)){
      if (y[i]==0) {
        yl[i] <- EnvStats::rnormTrunc(1, mean = X[i,]%*%bg, sd = 1, max = 0)
      }else if(y[i]==1){
        yl[i] <- EnvStats::rnormTrunc(1, mean = X[i,]%*%bg, sd = 1, min = 0, max = exp(tau[g,1]))
      }else if(y[i]==2){
        yl[i] <- EnvStats::rnormTrunc(1, mean = X[i,]%*%bg, sd = 1, min = exp(tau[g,1]), max = exp(tau[g,2])+exp(tau[g,1]))
      }else{
        yl[i] <- EnvStats::rnormTrunc(1, mean = X[i,]%*%bg, sd = 1, min = exp(tau[g,2])+exp(tau[g,1]))
      }
    }
    #Gibbs sampling for beta
    if(sum(is.nan(yl))>0 | sum(is.na(yl))>0 | sum(yl)==Inf){
      betasamples[g,] <- betasamples[g-1,]
    }else{
      b1 <- B1%*%(t(X)%*%yl + solve(B0)%*%b0)
      betasamples[g,] <- mvrnorm(1, mu = b1, Sigma = B1)
    }
    setWinProgressBar(pb, g, title=paste( round(g/G*100, 0),"% done"))
  }
  close(pb)
  return(cbind(betasamples, tau, ar))
}
###########################################Bayesian Estimation###################################
#Hyperparameters
d0 <- rep(0,2)
D0 <- diag(2)*10000
b0 <- rep(0,ncol(X))
B0 <- diag(ncol(X))*10000
#Estimation
param0 <- rep(0, ncol(X) + dplyr::n_distinct(y)-1)
G <- 1000
tun <- 1
VarProp <- tun*solve(solve(CovarML[4:5, 4:5]) + solve(D0))
param_sample <- mhop(param0, G)
#Burn in
B <- round(0.2*G)
param_sample <- param_sample[(B+1):G,]
mcmc0 <- coda::mcmc(param_sample[, 1:(ncol(X) + dplyr::n_distinct(y) - 2)])
summary(mcmc0)
summary(coda::mcmc(cbind(exp(param_sample[, 4]),exp(param_sample[, 4])+exp(param_sample[, 5]))))

########################## Negative binomial: Simulation ##########################
rm(list = ls())
set.seed(010101)
N <- 2000 # Sample size
x1 <- runif(N); x2 <- rnorm(N)
X <- cbind(1, x1, x2)
k <- dim(X)[2]
B <- rep(1, k)
alpha <- 1.2
gamma <- exp(alpha)
lambda <- exp(X%*%B)
y <- rnbinom(N, mu = lambda, size = gamma)
table(y)
# MCMC parameters
mcmc <- 10000
burnin <- 1000
thin <- 5
iter <- mcmc + burnin
keep <- seq(burnin, iter, thin)
sbeta <- 2.93/sqrt(k); salpha <- 2.93
# Hyperparameters: Priors
B0 <- 1000*diag(k); b0 <- rep(0, k)
alpha0 <- 0.5; delta0 <- 0.1
DataNB <- list(y = y, X = X)
mcmcNB <- list(R = mcmc, keep = thin, s_beta = sbeta, s_alpha = salpha)
PriorNB <- list(betabar = b0, A = solve(B0), a = alpha0, b = delta0)
ResultBayesm <- bayesm::rnegbinRw(Data = DataNB, Mcmc = mcmcNB, Prior = PriorNB)
summary(ResultBayesm$alphadraw)
summary(ResultBayesm$betadraw)

########################## Quantile regression: Value of soccer players ##########################
rm(list = ls())
set.seed(010101)
Data <- read.csv("DataApplications/1ValueFootballPlayers.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(Data)
y <- log(Value) 
X <- cbind(1, Perf, Age, Age2, NatTeam, Goals, Exp, Exp2)
RegLS <- lm(y ~ X -1)
k <- dim(X)[2]
N <- dim(X)[1]
# Hyperparameters
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
# MCMC parameters
mcmc <- 5000
burnin <- 1000
tot <- mcmc + burnin
thin <- 1
# Quantile
tau <- 0.5
theta <- (1-2*tau)/(tau*(1-tau))
psi2 <- 2/(tau*(1-tau))
an2 <- 2+theta^2/psi2
# Gibbs sampler
PostBeta <- function(e){
  Bn <- solve(B0i + psi2^(-1)*t(X)%*%diag(1/e)%*%X)
  bn <- Bn%*%(B0i%*%b0 + psi2^(-1)*t(X)%*%diag(1/e)%*%(y-theta*e))
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostE <- function(Beta, i){
  dn2 <-(y[i]-X[i,]%*%Beta)^2/psi2
  ei <- GIGrvg::rgig(1, chi = dn2, psi = an2, lambda = 1/2)
  return(ei)
}
PostBetas <- matrix(0, mcmc+burnin, k)
Beta <- RegLS$coefficients
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  e <- sapply(1:N, function(i){PostE(Beta = Beta, i)})
  Beta <- PostBeta(e = e)
  PostBetas[s,] <- Beta
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
colnames(PosteriorBetas) <- c("Intercept", "Perf", "Age", "Age2", "NatTeam", "Goals", "Exp", "Exp2")
summary(coda::mcmc(PosteriorBetas))

########################## Bayesian bootstrap: Simulation ##########################
rm(list = ls())
set.seed(010101)
N <- 1000 # Sample size
x1 <- runif(N); x2 <- rnorm(N)
X <- cbind(x1, x2)
k <- dim(X)[2]
B <- rep(1, k+1)
sig2 <- 1
u <- rnorm(N, 0, sig2)
y <- cbind(1, X)%*%B + u
data <- as.data.frame(cbind(y, X))
names(data) <- c("y", "x1", "x2")
Reg <- function(d){
  Reg <- lm(y ~ x1 + x2, data = d)
  Bhat <- Reg$coef
  return(Bhat)
}
Reg(data)
S <- 10000
BB <- bayesboot::bayesboot(data = data, statistic = Reg, R = S)
plot(BB)
