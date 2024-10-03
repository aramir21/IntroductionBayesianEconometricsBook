rm(list = ls())
set.seed(010101)
########################## Linear regression: Value of soccer players ##########################
Data <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/1ValueFootballPlayers.csv", sep = ",", header = TRUE, quote = "")
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
# MCMC parameters
mcmc <- 5000
burnin <- 5000
tot <- mcmc + burnin
thin <- 1
# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
posterior  <- MCMCpack::MCMCregress(y~X-1, b0=b0, B0 = B0i, c0 = a0, d0 = d0, burnin = burnin, mcmc = mcmc, thin = thin)
summary(coda::mcmc(posterior))
# Posterior distributions programming the Gibbs sampling
# Auxiliary parameters
XtX <- t(X)%*%X
bhat <- solve(XtX)%*%t(X)%*%y
an <- a0 + N
# Gibbs sampling functions
PostSig2 <- function(Beta){
  dn <- d0 + t(y - X%*%Beta)%*%(y - X%*%Beta)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBeta <- function(sig2){
  Bn <- solve(B0i + sig2^(-1)*XtX)
  bn <- Bn%*%(B0i%*%b0 + sig2^(-1)*XtX%*%bhat)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostBetas <- matrix(0, mcmc+burnin, k)
PostSigma2 <- rep(0, mcmc+burnin)
Beta <- rep(0, k)
for(s in 1:tot){
  sig2 <- PostSig2(Beta = Beta)
  PostSigma2[s] <- sig2
  Beta <- PostBeta(sig2 = sig2)
  PostBetas[s,] <- Beta
}
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
colnames(PosteriorBetas) <- c("Intercept", "Perf", "Age", "Age2", "NatTeam", "Goals", "Exp", "Exp2")
summary(coda::mcmc(PosteriorBetas))
PosteriorSigma2 <- PostSigma2[keep]
summary(coda::mcmc(PosteriorSigma2))
summary(coda::mcmc(exp(PosteriorBetas[,5])-1))

########################## Logit: Simulation ##########################
# Simulate data
rm(list = ls())
set.seed(010101)
N <- 10000 # Sample size
B <- c(0.5, 0.8, -1.2) # Population location parameters
x2 <- rnorm(N) # Regressor
x3 <- rnorm(N) # Regressor
X <- cbind(1, x2, x3) # Regressors
XB <- X%*%B
PY <- exp(XB)/(1 + exp(XB)) # Probability of Y = 1
Y <- rbinom(N, 1, PY) # Draw Y's
table(Y) # Frequency
# write.csv(cbind(Y, x2, x3), file = "DataSimulations/LogitSim.csv") # Export data

# MCMC parameters
iter <- 5000; burnin <- 1000; thin <- 5; tune <- 1
# Hyperparameters
K <- dim(X)[2] 
b0 <- rep(0, K)
c0 <- 1000
B0 <- c0*diag(K)
B0i <- solve(B0)
# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
RegLog <- MCMCpack::MCMClogit(Y~X-1, mcmc = iter, burnin = burnin, thin = thin, b0 = b0, B0 = B0i, tune = tune)
summary(RegLog)
# Posterior distributions programming the Metropolis-Hastings algorithm
# This function sets the M-H sampler using as default a hyperparameter mean equal to 0
# and a covariance equal to 1000 times a identity matrix, a tuning parameter equal to 1,
# 1000 post burn-in iterations, and the latter is equal to 500.
MHfunc <- function(y, X, b0 = rep(0, dim(X)[2] + 1), B0 = 1000*diag(dim(X)[2] + 1), tau = 1, 
                   iter = 6000, burnin = 1000, thin = 5){
  Xm <- cbind(1, X) # Regressors
  K <- dim(Xm)[2] # Number of location parameters
  BETAS <- matrix(0, iter + burnin, K) # Space for posterior chains
  Reg <- glm(y ~ Xm - 1, family = binomial(link = "logit")) # Maximum likelihood estimation
  BETA <- Reg$coefficients # Maximum likelihood parameter estimates 
  tot <- iter + burnin # Total iterations M-H algorithm
  COV <- vcov(Reg) # Maximum likelihood covariance matrix
  COVt <- tau^2*solve(solve(B0) + solve(COV)) # Covariance matrix for the proposal distribution
  Accep <- rep(0, tot) # Space for calculating the acceptance rate
  # create progress bar in case that you want to see iterations progress
  pb <- winProgressBar(title = "progress bar", min = 0,
                       max = tot, width = 300)
  for(it in 1:tot){
    BETAc <- BETA + MASS::mvrnorm(n = 1, mu = rep(0, K), Sigma = COVt) # Candidate location parameter
    likecand <- sum((Xm%*%BETAc) * Y - apply(Xm%*%BETAc, 1, function(x) log(1 + exp(x)))) # Log likelihood for the candidate
    likepast <- sum((Xm%*%BETA) * Y - apply((Xm%*%BETA), 1, function(x) log(1 + exp(x)))) # Log likelihood for the actual draw
    priorcand <- (-1/2)*crossprod((BETAc - b0), solve(B0))%*%(BETAc - b0) # Log prior for candidate
    priorpast <- (-1/2)*crossprod((BETA - b0), solve(B0))%*%(BETA - b0) # Log prior for actual draw
    alpha <- min(1, exp((likecand + priorcand) - (likepast + priorpast))) #Probability of selecting candidate
    u <- runif(1) # Decision rule for selecting candidate
    if(u < alpha){
      BETA <- BETAc # Changing reference for candidate if selected
      Accep[it] <- 1 # Indicator if the candidate is accepted
    } 
    BETAS[it, ] <- BETA # Saving draws
    setWinProgressBar(pb, it, title=paste( round(it/tot*100, 0),
                                          "% done"))
  }
  close(pb)
  keep <- seq(burnin, tot, thin)
  return(list(Bs = BETAS[keep[-1], ], AceptRate = mean(Accep[keep[-1]])))
}
Posterior <- MHfunc(y = Y, X = cbind(x2, x3), iter = iter, burnin = burnin, thin = thin) # Running our M-H function changing some default parameters.
paste("Acceptance rate equal to", round(Posterior$AceptRate, 2), sep = " ")
PostPar <- coda::mcmc(Posterior$Bs)
# Names
colnames(PostPar) <- c("Cte", "x1", "x2")
# Summary posterior draws
summary(PostPar)
# Trace and density plots
plot(PostPar)
# Autocorrelation plots
coda::autocorr.plot(PostPar)
# Convergence diagnostics
coda::geweke.diag(PostPar)
coda::raftery.diag(PostPar,q=0.5,r=0.05,s = 0.95)
coda::heidel.diag(PostPar)

########################## Probit: Hospitalization ##########################
rm(list = ls())
set.seed(010101)
mydata <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/2HealthMed.csv", sep = ",", header = TRUE, quote = "")
attach(mydata)
str(mydata)
y <- Hosp # Dependent variables
X <- cbind(1, SHI, Female, Age, Age2, Est2, Est3, Fair, Good, Excellent) # Regressors
K <- 10 # Number of regressors
# Hyperparameters
b0 <- rep(0, K) # Prio mean
B0i <- diag(K) # Prior precision (inverse of covariance)
Prior <- list(betabar = b0, A = B0i) # Prior list
Data <- list(y = y, X = X) # Data list
Mcmc <- list(R = 10000, keep = 1, nprint = 0) # MCMC parameters
RegProb <- bayesm::rbprobitGibbs(Data = Data, Prior = Prior, Mcmc = Mcmc) # Inference using bayesm package
PostPar <- coda::mcmc(RegProb$betadraw) # Posterior draws
colnames(PostPar) <- c("Cte", "SHI", "Female", "Age", "Age2", "Est2",
                       "Est3", "Fair", "Good", "Excellent") # Names
summary(PostPar) # Posterior summary

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
y1<-y==1
y2<-y==2
y3<-y==3
# Log likelihood
log.L<- function(Beta){
  V1<-Beta[1]+Beta[3]*X4+X1%*%Beta[5:7]
  V2<-Beta[2]+Beta[4]*X4+X2%*%Beta[5:7]
  V3<- X3%*%Beta[5:7]
  suma<-exp(V1)+exp(V2)+exp(V3)
  p11<-exp(V1)/suma
  p22<-exp(V2)/suma
  p33<-exp(V3)/suma
  suma2<-NULL
  for(i in 1:N){
    suma1<-y1[i]*log(p11[i])+y2[i]*log(p22[i])+y3[i]*log(p33[i])
    suma2<-c(suma2,suma1)}
  logL<-sum(suma2)
  return(-logL)
}
# Parameters: Proposal
k <- 7
res.optim<-optim(rep(0, k), log.L, method="BFGS", hessian=TRUE)
MeanT <- res.optim$par
ScaleT <- as.matrix(Matrix::forceSymmetric(solve(res.optim$hessian))) # Force this matrix to be symmetric
# Hyperparameters: Priors
B0 <- 1000*diag(k)
b0 <- rep(0, k)

MHfunction <- function(iter, tuning){
  Beta <- rep(0, k)
  Acept <- NULL
  BetasPost <- matrix(NA, iter, k)
  pb <- winProgressBar(title = "progress bar", min = 0, max = iter, width = 300)
  for(s in 1:iter){
    LogPostBeta <- -log.L(Beta) + mvtnorm::dmvnorm(Beta, mean = b0, sigma = B0, log = TRUE)
    BetaC <- c(LaplacesDemon::rmvt(n=1, mu = MeanT, S = ScaleT, df = tuning))
    LogPostBetaC <- -log.L(BetaC) + mvtnorm::dmvnorm(BetaC, mean = b0, sigma = B0, log = TRUE)
    alpha <- min(exp((LogPostBetaC-mvtnorm::dmvt(BetaC, delta = MeanT, sigma = ScaleT, df = tuning, log = TRUE))-(LogPostBeta-mvtnorm::dmvt(Beta, delta = MeanT, sigma = ScaleT, df = tuning, log = TRUE))) ,1)
    u <- runif(1)
    if(u <= alpha){
      Acepti <- 1
      Beta <- BetaC
    }else{
      Acepti <- 0
      Beta <- Beta
    }
    BetasPost[s, ] <- Beta
    Acept <- c(Acept, Acepti)
    setWinProgressBar(pb, s, title=paste( round(s/iter*100, 0),"% done"))
  }
  close(pb)
  AcepRate <- mean(Acept)
  Results <- list(AcepRate = AcepRate, BetasPost = BetasPost)
  return(Results)
}
# MCMC parameters
mcmc <- 10000
burnin <- 1000
thin <- 5
iter <- mcmc + burnin
keep <- seq(burnin, iter, thin)
tuning <- 6 # Degrees of freedom of multivariate Student's t
# Run M-H
ResultsPost <- MHfunction(iter = iter, tuning = tuning)
ResultsPost$AcepRate
summary(coda::mcmc(ResultsPost$BetasPost[keep[-1], ]))

########################## Ordered probit: Preventive doctor visits ##########################
rm(list = ls())
set.seed(010101)
Data <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/2HealthMed.csv", sep = ",", header = TRUE, quote = "")
attach(Data)
y <- MedVisPrevOr 
# MedVisPrevOr: Oredered variable for preventive visits to doctors in one year: 1 (none), 2 (once), ... 6 (five or more)
X <- cbind(SHI, Female, Age, Age2, Est2, Est3, Fair, Good, Excellent, PriEd, HighEd, VocEd, UnivEd)
# Subsideized health care program: SHI
# Sex: Female
# Age: Linear and squared
# Socioeconomic condition: Est2 and Est3, the lowest (Est1) is the baseline category
# Health perception: Fair, Good, Excellent 
# Education levels: PriEd, HighEd, VocEd, UnivEd
k <- dim(X)[2]
L <- length(table(y))
# Hyperparameters
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
gamma0 <- rep(0, L-2)
Gamma0 <- diag(L-2)
# MCMC parameters
mcmc <- 60000+1
thin <- 5
tuningPar <- 1/(L-2)^0.5
DataApp <- list(y = y, X = X, k = L)
Prior <- list(betabar = b0, A = solve(B0), dstarbar = gamma0, Ad = Gamma0)
mcmcpar <- list(R = mcmc, keep = 5, s = tuningPar)
PostBeta <- bayesm::rordprobitGibbs(Data = DataApp, Prior = Prior, Mcmc = mcmcpar)
BetasPost <- coda::mcmc(PostBeta[["betadraw"]])
colnames(BetasPost) <- c("SHI", "Female", "Age", "Age2", "Est2", "Est3", "Fair", "Good", "Excellent", "PriEd", "HighEd", "VocEd", "UnivEd")
summary(BetasPost)
# Convergence diagnostics
coda::geweke.diag(BetasPost)
coda::raftery.diag(BetasPost,q=0.5,r=0.05,s = 0.95)
coda::heidel.diag(BetasPost)
# Cut offs
Cutoffs <- PostBeta[["cutdraw"]]
summary(Cutoffs)
coda::geweke.diag(Cutoffs)
coda::heidel.diag(Cutoffs)
coda::raftery.diag(Cutoffs[,-1],q=0.5,r=0.05,s = 0.95)

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
# log likelihood
logLik <- function(par){
  alpha <- par[1]; beta <- par[2:(k+1)]
  gamma <- exp(alpha)
  lambda <- exp(X%*%beta)
  logLikNB <- sum(sapply(1:N, function(i){dnbinom(y[i], size = gamma, mu = lambda[i], log = TRUE)}))
  return(-logLikNB)
}
# Parameters: Proposal
par0 <- rep(0.5, k+1)
res.optim <- optim(par0, logLik, method="BFGS", hessian=TRUE)
res.optim$par
res.optim$convergence
Covar <- solve(res.optim$hessian) 
CovarBetas <- Covar[2:(k+1),2:(k+1)]
VarAlpha <- Covar[1:1]
# Hyperparameters: Priors
B0 <- 1000*diag(k); b0 <- rep(0, k)
alpha0 <- 0.5; delta0 <- 0.1

MHfunction <- function(iter, sbeta, salpha){
  Beta <- rep(0, k)
  Acept1 <- NULL
  Acept2 <- NULL
  BetasPost <- matrix(NA, iter, k)
  alpha <- 1
  alphaPost <- rep(NA, iter)
  par <- c(alpha, Beta)
  pb <- winProgressBar(title = "progress bar", min = 0, max = iter, width = 300)
  for(s in 1:iter){
    LogPostBeta <- -logLik(par) + dgamma(alpha, shape = alpha0, scale = delta0, log = TRUE) + mvtnorm::dmvnorm(Beta, mean = b0, sigma = B0, log = TRUE)
    BetaC <- c(MASS::mvrnorm(1, mu = Beta, Sigma = sbeta^2*CovarBetas))
    parC <- c(alpha, BetaC)
    LogPostBetaC <- -logLik(parC) + dgamma(alpha, shape = alpha0, scale = delta0, log = TRUE) +  mvtnorm::dmvnorm(BetaC, mean = b0, sigma = B0, log = TRUE)
    alpha1 <- min(exp((LogPostBetaC - mvtnorm::dmvnorm(BetaC, mean = Beta, sigma = sbeta^2*CovarBetas, log = TRUE))-(LogPostBeta - mvtnorm::dmvnorm(Beta, mean = Beta, sigma = sbeta^2*CovarBetas, log = TRUE))),1)
    u1 <- runif(1)
    if(u1 <= alpha1){
      Acept1i <- 1
      Beta <- BetaC
    }else{
      Acept1i <- 0
      Beta <- Beta
    }
    par <- c(alpha, Beta)
    LogPostBeta <- -logLik(par) + dgamma(alpha, shape = alpha0, scale = delta0, log = TRUE) + mvtnorm::dmvnorm(Beta, mean = b0, sigma = B0, log = TRUE)
    alphaC <- rnorm(1, mean = alpha, sd = salpha*VarAlpha^0.5)
    parC <- c(alphaC, Beta)
    LogPostBetaC <- -logLik(parC) + dgamma(alphaC, shape = alpha0, scale = delta0, log = TRUE) +  mvtnorm::dmvnorm(Beta, mean = b0, sigma = B0, log = TRUE)
    alpha2 <- min(exp((LogPostBetaC - dnorm(alphaC, mean = alpha, sd = salpha*VarAlpha^0.5, log = TRUE))-(LogPostBeta - dnorm(alpha, mean = alpha, sd = salpha*VarAlpha^0.5, log = TRUE))),1)
    u2 <- runif(1)
    if(u2 <= alpha2){
      Acept2i <- 1
      alpha <- alphaC
    }else{
      Acept2i <- 0
      alpha <- alpha
    }
    
    BetasPost[s, ] <- Beta
    alphaPost[s] <- alpha
    Acept1 <- c(Acept1, Acept1i)
    Acept2 <- c(Acept2, Acept2i)
    setWinProgressBar(pb, s, title=paste( round(s/iter*100, 0),"% done"))
  }
  close(pb)
  AcepRateBeta <- mean(Acept1); AcepRateAlpha <- mean(Acept2)
  Results <- list(AcepRateBeta = AcepRateBeta, AcepRateAlpha = AcepRateAlpha, BetasPost = BetasPost, alphaPost = alphaPost)
  return(Results)
}
# MCMC parameters
mcmc <- 10000
burnin <- 1000
thin <- 5
iter <- mcmc + burnin
keep <- seq(burnin, iter, thin)
sbeta <- 2.93/sqrt(k); salpha <- 2.93
# Run M-H
ResultsPost <- MHfunction(iter = iter, sbeta = sbeta, salpha = salpha)
ResultsPost$AcepRateBeta
ResultsPost$AcepRateAlpha
summary(coda::mcmc(ResultsPost$BetasPost[keep[-1], ]))
summary(coda::mcmc(ResultsPost$alphaPost[keep[-1]]))

DataNB <- list(y = y, X = X)
mcmcNB <- list(R = mcmc, keep = thin, s_beta = sbeta, s_alpha = salpha)
PriorNB <- list(betabar = b0, A = solve(B0), a = alpha0, b = delta0)
ResultBayesm <- bayesm::rnegbinRw(Data = DataNB, Mcmc = mcmcNB, Prior = PriorNB)
summary(ResultBayesm$alphadraw)
summary(ResultBayesm$betadraw)

########################## Tobit model: Value of soccer players ##########################
rm(list = ls())
set.seed(010101)
Data <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/1ValueFootballPlayers.csv", sep = ",", header = TRUE, quote = "")
attach(Data)
y <- log(ValueCens) 
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
an <- a0 + N
# MCMC parameters
mcmc <- 50000
burnin <- 10000
tot <- mcmc + burnin
thin <- 1
# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
posterior  <- MCMCpack::MCMCtobit(y~X-1, b0=b0, B0 = B0i, c0 = a0, d0 = d0, burnin = burnin, mcmc = mcmc, thin = thin, below = 13.82, above = Inf)
summary(coda::mcmc(posterior))
# Posterior distributions programming the Gibbs sampling
# Auxiliary parameters
XtX <- t(X)%*%X
# Gibbs sampling functions
PostBeta <- function(Yl, sig2){
  Bn <- solve(B0i + sig2^(-1)*XtX)
  bn <- Bn%*%(B0i%*%b0 + sig2^(-1)*t(X)%*%Yl)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostYl <- function(Beta, L, U, i){
  Ylmean <- X[i,]%*%Beta
  if(y[i] == L){
    Yli <- truncnorm::rtruncnorm(1, a = -Inf, b = L, mean = Ylmean, sd = sig2^0.5)
  }else{
    if(y[i] == U){
      Yli <- truncnorm::rtruncnorm(1, a = U, b = Inf, mean = Ylmean, sd = sig2^0.5)
    }else{
      Yli <- y[i]
    }
  }
  return(Yli)
}
PostSig2 <- function(Beta, Yl){
  dn <- d0 + t(Yl - X%*%Beta)%*%(Yl - X%*%Beta)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}
PostBetas <- matrix(0, mcmc+burnin, k)
Beta <- rep(0, k)
PostSigma2 <- rep(0, mcmc+burnin)
sig2 <- 1
L <- log(1000000)
U <- Inf
# create progress bar in case that you want to see iterations progress
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Yl <- sapply(1:N, function(i){PostYl(Beta = Beta, L = L, U = U, i)})
  Beta <- PostBeta(Yl = Yl, sig2)
  sig2 <- PostSig2(Beta = Beta, Yl = Yl) 
  PostBetas[s,] <- Beta
  PostSigma2[s] <- sig2
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0), "% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
colnames(PosteriorBetas) <- c("Intercept", "Perf", "Age", "Age2", "NatTeam", "Goals", "Exp", "Exp2")
summary(coda::mcmc(PosteriorBetas))
summary(coda::mcmc(PostSigma2[keep]))

######################### Quantile model: Value of soccer players ##########################
rm(list = ls())
set.seed(010101)
Data <- read.csv("DataApplications/1ValueFootballPlayers.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(Data)
y <- log(ValueCens) 
# Value: Market value in Euros (2017) of soccer players
# Regressors quantity including intercept
X <- cbind(1, Perf, Age, Age2, NatTeam, Goals, Exp, Exp2)
# Perf: Performance. Perf2: Performance squared. Age: Age; Age: Age squared. 
# NatTeam: Indicator of national team. Goals: Scored goals. Goals2: Scored goals squared
# Exp: Years of experience. Exp2: Years of experience squared. Assists: Number of assists
k <- dim(X)[2]
N <- dim(X)[1]
# Hyperparameters
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
# MCMC parameters
mcmc <- 50000
burnin <- 10000
tot <- mcmc + burnin
thin <- 1
# Quantile
q <- 0.5
# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
posterior  <- MCMCpack::MCMCquantreg(y~X-1, tau = q, b0=b0, B0 = B0i, burnin = burnin, mcmc = mcmc, thin = thin, below = 13.82, above = Inf)
summary(coda::mcmc(posterior))
q <- 0.9
# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
posterior  <- MCMCpack::MCMCquantreg(y~X-1, tau = q, b0=b0, B0 = B0i, burnin = burnin, mcmc = mcmc, thin = thin, below = 13.82, above = Inf)
summary(coda::mcmc(posterior))

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
S <- 10000; alpha <- 1
BB <- function(S, df, alpha){
  Betas <- matrix(NA, S, dim(df)[2])
  N <- dim(df)[1]
  pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
  for(s in 1:S){
    g <- LaplacesDemon::rdirichlet(N, alpha)
    ids <- sample(1:N, size = N, replace = TRUE, prob = g)
    datas <- df[ids,]
    names(datas) <- names(df)
    Bs <- Reg(d = datas)
    Betas[s, ] <- Bs
    setWinProgressBar(pb, s, title=paste( round(s/S*100, 0), "% done"))
  }
  close(pb)
  return(Betas)
}
BBs <- BB(S = S, df = data, alpha = alpha)
summary(coda::mcmc(BBs))
