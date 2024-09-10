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
d0 <- 0.001/2
a0 <- 0.001/2
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
  taui <- rgamma(1, v1n, v2n)
  return(taui)
}
PostBetas <- matrix(0, mcmc+burnin, k)
PostSigma2 <- rep(0, mcmc+burnin)
Beta <- rep(0, k)
tau <- rep(1, N)
for(s in 1:tot){
  sig2 <- PostSig2(Beta = Beta, tau = tau)
  PostSigma2[s] <- sig2
  Beta <- PostBeta(sig2 = sig2, tau = tau)
  PostBetas[s,] <- Beta
  tau <- sapply(1:N, function(i){PostTau(sig2 = sig2, Beta = Beta, i)})
}
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
dat <-data.frame(mode,X1[,1],X2[,1],X3[,1],X1[,2],X2[,2],X3[,2],X1[,3],X2[,3],X3[,3],X4)
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

########################## Ordered probit application: Preventive visits to physicians ##########################


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
