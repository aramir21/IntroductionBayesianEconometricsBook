########################## Effects of capital on productivity USA states: Longitudinal normal model ########################## 
rm(list = ls())
set.seed(12345)
DataGSP <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/8PublicCap.csv", sep = ",", header = TRUE, quote = "")
attach(DataGSP)
K1 <- 5; K2 <- 1
b0 <- rep(0, K1); B0 <- diag(K1)
r0 <- 5; R0 <- diag(K2)
a0 <- 0.001; d0 <- 0.001
Resultshreg <- MCMCpack::MCMChregress(fixed = log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, random = ~1, 
                                      group = "id", data = DataGSP, burnin = 5000, mcmc = 10000, thin = 1,
                                      r = r0, R = R0, nu = a0, delta = d0)
Betas <- Resultshreg[["mcmc"]][,1:K1]
Sigma2RanEff <- Resultshreg[["mcmc"]][,54]
Sigma2 <- Resultshreg[["mcmc"]][,55]
summary(Betas)
summary(Sigma2RanEff)
summary(Sigma2)
summary(Sigma2RanEff/(Sigma2RanEff+Sigma2))
# Convergence diagnostics
coda::geweke.diag(Betas)
coda::raftery.diag(Betas,q=0.5,r=0.05,s = 0.95)
coda::heidel.diag(Betas)

########################## Simulation exercise: Longitudinal normal model ########################## 
rm(list = ls())
set.seed(010101)
NT <- 2000
N <- 50
id <- c(1:N, sample(1:N, NT - N,replace=TRUE))
table(id)
x1 <- rnorm(NT); x2 <- rnorm(NT); x3 <- rnorm(NT) 
X <- cbind(1, x1, x2, x3)
K1 <- dim(X)[2]
w1 <- rnorm(NT) 
W <- cbind(1, w1)
K2 <- dim(W)[2]
B <- c(0.5, 0.4, 0.6, -0.6)
D <- c(0.7, 0.6)
b1 <- rnorm(N, 0, sd = D[1]^0.5)
b2 <- rnorm(N, 0, sd = D[2]^0.5)
b <- cbind(b1, b2)
sig2 <- 0.1
u <- rnorm(NT, 0, sd = sig2^0.5)
y <- NULL
for(i in 1:NT){
  yi <- X[i,]%*%B + W[i,]%*%b[id[i],] + u[i] 
  y <- c(y, yi)
}
Data <- as.data.frame(cbind(y, x1, x2, x3, w1, id))
mcmc <- 5000; burnin <- 1000; thin <- 1; tot <- mcmc + burnin
b0 <- rep(0, K1); B0 <- diag(K1); B0i <- solve(B0) 
r0 <- K2; R0 <- diag(K2); a0 <- 0.001; d0 <- 0.001
# MCMChregress
Resultshreg <- MCMCpack::MCMChregress(fixed = y~x1 + x2 + x3, random = ~w1, group="id",
                      data = Data, burnin = burnin, mcmc = mcmc, thin = thin, 
                      mubeta = b0, Vbeta = B0,
                      r = r0, R = R0, nu = a0, delta = d0)
Betas <- Resultshreg[["mcmc"]][,1:K1]
Sigma2RanEff <- Resultshreg[["mcmc"]][,c(K2*N+K1+1, 2*N+K1+K2^2)]
Sigma2 <- Resultshreg[["mcmc"]][,K2*N+K1+K2^2+1]
summary(Betas)
summary(Sigma2RanEff)
summary(Sigma2)
coda::geweke.diag(Betas)
coda::raftery.diag(Betas,q=0.5, r=0.05, s = 0.95)
coda::heidel.diag(Betas)
# Gibbs by hand
PostBeta <- function(sig2, D){
  XVX <- matrix(0, K1, K1)
  XVy <- matrix(0, K1, 1)
  for(i in 1:N){
    ids <- which(id == i)
    Ti <- length(ids)
    Wi <- W[ids, ]
    Vi <- sig2*diag(Ti) + Wi%*%D%*%t(Wi)
    ViInv <- solve(Vi)
    Xi <- X[ids, ]
    XVXi <- t(Xi)%*%ViInv%*%Xi
    XVX <- XVX + XVXi
    yi <- y[ids]
    XVyi <- t(Xi)%*%ViInv%*%yi
    XVy <- XVy + XVyi
  }
  Bn <- solve(B0i + XVX)
  bn <- Bn%*%(B0i%*%b0 + XVy)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
# This version of the function follows MCMCpack function named cMCMChregress.cc that uses 
# http://en.wikipedia.org/wiki/Woodbury_matrix_identity with A=(1/V_run)*Id
PostBetaNew <- function(sig2, D){
  Di <- solve(D)
  XVXn <- matrix(0, K1, K1)
  XVyn <- matrix(0, K1, 1)
  for(i in 1:N){
    ids <- which(id == i)
    Ti <- length(ids)
    Wi <- W[ids, ]
    Xi <- X[ids, ]
    yi <- y[ids]
    sig2i <- sig2^(-1)
    XVXi <- t(Xi)%*%Xi*sig2i - sig2i^(2)*t(Xi)%*%Wi%*%solve(Di+sig2i*t(Wi)%*%Wi)%*%t(Wi)%*%Xi
    XVyi <- t(Xi)%*%yi*sig2i - sig2i^(2)*t(Xi)%*%Wi%*%solve(Di+sig2i*t(Wi)%*%Wi)%*%t(Wi)%*%yi 
    XVXn <- XVXn + XVXi
    XVyn <- XVyn + XVyi
  }
  BnNew <- solve(B0i + XVXn)
  bnNew <- BnNew%*%(B0i%*%b0 + XVyn)
  Beta <- MASS::mvrnorm(1, bnNew, BnNew)
  return(Beta)
}
Postb <- function(Beta, sig2, D){
  Di <- solve(D)
  bis <- matrix(0, N, K2)
  for(i in 1:N){
    ids <- which(id == i)
    Wi <- W[ids, ]
    Xi <- X[ids, ]
    yi <- y[ids]
    Wtei <- sig2^(-1)*t(Wi)%*%(yi - Xi%*%Beta)
    Bni <- solve(sig2^(-1)*t(Wi)%*%Wi + Di)
    bni <- Bni%*%Wtei
    bi <- MASS::mvrnorm(1, bni, Bni)
    bis[i, ] <- bi
  }
  return(bis)
}
PostSig2 <- function(Beta, bs){
  an <- a0 + 0.5*NT
  ete <- 0
  for(i in 1:N){
    ids <- which(id == i)
    Xi <- X[ids, ]
    yi <- y[ids]
    Wi <- W[ids, ]
    ei <- yi - Xi%*%Beta - Wi%*%bs[i, ]
    etei <- t(ei)%*%ei
    ete <- ete + etei
  }
  dn <- d0 + 0.5*ete 
  sig2 <- MCMCpack::rinvgamma(1, shape = an, scale = dn)
  return(sig2)
}
PostD <- function(bs){
  rn <- r0 + N
  btb <- matrix(0, K2, K2)
  for(i in 1:N){
    bsi <- bs[i, ]
    btbi <- bsi%*%t(bsi)
    btb <- btb + btbi
  }
  Rn <- d0*R0 + btb
  Sigma <- MCMCpack::riwish(v = rn, S = Rn)
  return(Sigma)
}
PostBetas <- matrix(0, tot, K1)
PostDs <- matrix(0, tot, K2*(K2+1)/2)
PostSig2s <- rep(0, tot)
Postbs <- array(0, c(N, K2, tot))
RegLS <- lm(y ~ X - 1)
SumLS <- summary(RegLS)
Beta <- SumLS[["coefficients"]][,1]
sig2 <- SumLS[["sigma"]]^2
D <- diag(K2)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  bs <- Postb(Beta = Beta, sig2 = sig2, D = D)
  D <- PostD(bs = bs)
  Beta <- PostBeta(sig2 = sig2, D = D)
  # Beta <- PostBetaNew(sig2 = sig2, D = D)
  sig2 <- PostSig2(Beta = Beta, bs = bs)
  PostBetas[s,] <- Beta
  PostDs[s,] <- matrixcalc::vech(D)
  PostSig2s[s] <- sig2
  Postbs[, , s] <- bs
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
Bs <- PostBetas[keep,]
Ds <- PostDs[keep,]
bs <- Postbs[, , keep]
sig2s <- PostSig2s[keep]
summary(coda::mcmc(Bs))
plot(coda::mcmc(Bs))
summary(coda::mcmc(Ds))
plot(coda::mcmc(Ds))
summary(coda::mcmc(sig2s))
plot(coda::mcmc(sig2s))
# Convergence diagnostics
coda::geweke.diag(Bs)
coda::raftery.diag(Bs, q = 0.5, r = 0.05, s = 0.95)
coda::heidel.diag(Bs)

########################## Simulation exercise: Longitudinal normal model with heteroskedasticity ########################## 
rm(list = ls())
set.seed(010101)
NT <- 2000
N <- 50
id <- c(1:N, sample(1:N, NT - N,replace=TRUE))
table(id)
x1 <- rnorm(NT); x2 <- rnorm(NT); x3 <- rnorm(NT) 
X <- cbind(1, x1, x2, x3)
K1 <- dim(X)[2]
w1 <- rnorm(NT) 
W <- cbind(1, w1)
K2 <- dim(W)[2]
B <- c(0.5, 0.4, 0.6, -0.6)
D <- c(0.7, 0.6)
b1 <- rnorm(N, 0, sd = D[1]^0.5)
b2 <- rnorm(N, 0, sd = D[2]^0.5)
b <- cbind(b1, b2)
v <- 5
tau <- rgamma(NT, shape = v/2, rate = v/2)
sig2 <- 0.1
u <- rnorm(NT, 0, sd = (sig2/tau)^0.5)
y <- NULL
for(i in 1:NT){
  yi <- X[i,]%*%B + W[i,]%*%b[id[i],] + u[i] 
  y <- c(y, yi)
}
Data <- as.data.frame(cbind(y, x1, x2, x3, w1, id))
mcmc <- 5000; burnin <- 1000; thin <- 1; tot <- mcmc + burnin
b0 <- rep(0, K1); B0 <- diag(K1); B0i <- solve(B0) 
r0 <- K2; R0 <- diag(K2); a0 <- 0.001; d0 <- 0.001
# MCMChregress
Resultshreg <- MCMCpack::MCMChregress(fixed = y~x1 + x2 + x3, random = ~w1, group="id",
                                      data = Data, burnin = burnin, mcmc = mcmc, thin = thin, 
                                      mubeta = b0, Vbeta = B0,
                                      r = r0, R = R0, nu = a0, delta = d0)
Betas <- Resultshreg[["mcmc"]][,1:K1]
Sigma2RanEff <- Resultshreg[["mcmc"]][,c(K2*N+K1+1, 2*N+K1+K2^2)]
Sigma2 <- Resultshreg[["mcmc"]][,K2*N+K1+K2^2+1]
summary(Betas)
summary(Sigma2RanEff)
summary(Sigma2)
# Gibbs by hand
PostBeta <- function(sig2, D, tau){
  XVX <- matrix(0, K1, K1)
  XVy <- matrix(0, K1, 1)
  for(i in 1:N){
    ids <- which(id == i)
    Ti <- length(ids)
    Wi <- W[ids, ]
    taui <- tau[ids]
    Vi <- sig2*solve(diag(1/taui)) + Wi%*%D%*%t(Wi)
    ViInv <- solve(Vi)
    Xi <- X[ids, ]
    XVXi <- t(Xi)%*%ViInv%*%Xi
    XVX <- XVX + XVXi
    yi <- y[ids]
    XVyi <- t(Xi)%*%ViInv%*%yi
    XVy <- XVy + XVyi
  }
  Bn <- solve(B0i + XVX)
  bn <- Bn%*%(B0i%*%b0 + XVy)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
# sig2 <- sig2; D <- matrix(c(0.7,0,0,0.6),2,2); tau <- tau
# rowMeans(replicate(100, PostBeta(sig2, D, tau)))
Postb <- function(Beta, sig2, D, tau){
  Di <- solve(D)
  bis <- matrix(0, N, K2)
  for(i in 1:N){
    ids <- which(id == i)
    Wi <- W[ids, ]
    Xi <- X[ids, ]
    yi <- y[ids]
    taui <- tau[ids]
    Taui <- solve(diag(1/taui))
    Wtei <- sig2^(-1)*t(Wi)%*%Taui%*%(yi - Xi%*%Beta)
    Bni <- solve(sig2^(-1)*t(Wi)%*%Taui%*%Wi + Di)
    bni <- Bni%*%Wtei
    bi <- MASS::mvrnorm(1, bni, Bni)
    bis[i, ] <- bi
  }
  return(bis)
}
# sig2 <- sig2; D <- matrix(c(0.7,0,0,0.6),2,2); tau <- tau; Beta <- B
# Ens <- replicate(100, rowMeans(Postb(Beta, sig2, D, tau) - b))
PostSig2 <- function(Beta, bs, tau){
  an <- a0 + 0.5*NT
  ete <- 0
  for(i in 1:N){
    ids <- which(id == i)
    Xi <- X[ids, ]
    yi <- y[ids]
    Wi <- W[ids, ]
    taui <- tau[ids]
    Taui <- solve(diag(1/taui))
    ei <- yi - Xi%*%Beta - Wi%*%bs[i, ]
    etei <- t(ei)%*%Taui%*%ei
    ete <- ete + etei
  }
  dn <- d0 + 0.5*ete 
  sig2 <- MCMCpack::rinvgamma(1, shape = an, scale = dn)
  return(sig2)
}
# tau <- tau; Beta <- B; bs <- b
# mean(replicate(100, PostSig2(Beta, bs, tau)))
PostD <- function(bs){
  rn <- r0 + N
  btb <- matrix(0, K2, K2)
  for(i in 1:N){
    bsi <- bs[i, ]
    btbi <- bsi%*%t(bsi)
    btb <- btb + btbi
  }
  Rn <- d0*R0 + btb
  Sigma <- MCMCpack::riwish(v = rn, S = Rn)
  return(Sigma)
}
# bs <- b
# rowMeans(replicate(100, PostD(bs)[c(1,4)]))
PostTau <- function(sig2, Beta, bs){
  v1n <- v + 1
  v2n <- NULL
  for(i in 1:NT){
    Xi <- X[i, ]
    yi <- y[i]
    Wi <- W[i, ]
    bi <- bs[id[i],]
    v2ni <- v + sig2^(-1)*(yi - Xi%*%Beta - Wi%*%bi)^2
    v2n <- c(v2n, v2ni)
  }
  tau <- rgamma(NT, shape = rep(v1n/2, NT), rate = v2n/2)
  return(tau)
}
# bs <- b; sig2 <- sig2; Beta <- B
# Enstau <- rowMeans(replicate(100, PostTau(sig2, Beta, bs)))
# summary(Enstau)
# summary(tau)
PostBetas <- matrix(0, tot, K1)
PostDs <- matrix(0, tot, K2*(K2+1)/2)
PostSig2s <- rep(0, tot)
Postbs <- array(0, c(N, K2, tot))
PostTaus <- matrix(0, tot, NT) 
RegLS <- lm(y ~ X - 1)
SumLS <- summary(RegLS)
Beta <- SumLS[["coefficients"]][,1]
sig2 <- SumLS[["sigma"]]^2
D <- diag(K2)
tau <- rgamma(NT, shape = v/2, rate = v/2) 
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  bs <- Postb(Beta = Beta, sig2 = sig2, D = D, tau = tau)
  D <- PostD(bs = bs)
  Beta <- PostBeta(sig2 = sig2, D = D, tau = tau)
  sig2 <- PostSig2(Beta = Beta, bs = bs, tau = tau)
  tau <- PostTau(sig2 = sig2, Beta = Beta, bs = bs)
  PostBetas[s,] <- Beta
  PostDs[s,] <- matrixcalc::vech(D)
  PostSig2s[s] <- sig2
  Postbs[, , s] <- bs
  PostTaus[s,] <- tau
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
Bs <- PostBetas[keep,]
Ds <- PostDs[keep,]
bs <- Postbs[, , keep]
sig2s <- PostSig2s[keep]
taus <- PostTaus[keep,]
summary(coda::mcmc(Bs))
plot(coda::mcmc(Bs))
summary(coda::mcmc(Ds))
plot(coda::mcmc(Ds))
summary(coda::mcmc(sig2s))
plot(coda::mcmc(sig2s))
summary(coda::mcmc(taus))
# Convergence diagnostics
coda::geweke.diag(Bs)
coda::raftery.diag(Bs, q = 0.5, r = 0.05, s = 0.95)
coda::heidel.diag(Bs)

########################## Visits to doctor: Longitudinal logit model ########################## 
rm(list = ls())
set.seed(12345)
Data <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/9VisitDoc.csv", sep = ",", header = TRUE, quote = "")
attach(Data)
K1 <- 7; K2 <- 2; N <- 9197
b0 <- rep(0, K1); B0 <- diag(K1)
r0 <- 5; R0 <- diag(K2)
a0 <- 0.001; d0 <- 0.001
RegLogit <- glm(DocVis ~ Age + Male + Sport + LogInc + GoodHealth + BadHealth, family = binomial(link = "logit"))
SumLogit <- summary(RegLogit)
Beta0 <- SumLogit[["coefficients"]][,1]
mcmc <- 10000; burnin <- 1000; thin <- 10
# MCMChlogit
Resultshlogit <- MCMCpack::MCMChlogit(fixed = DocVis ~ Age + Male + Sport + LogInc + GoodHealth + BadHealth, random = ~Sozh, group="id",
                                      data = Data, burnin = burnin, mcmc = mcmc, thin = thin, 
                                      mubeta = b0, Vbeta = B0,
                                      r = r0, R = R0, nu = a0, delta = d0,
                                      beta.start = Beta0, FixOD = 1)

Betas <- Resultshlogit[["mcmc"]][,1:K1]
Sigma2RanEff <- Resultshlogit[["mcmc"]][,c(K2*N+K1+1, 2*N+K1+K2^2)]
summary(Betas)
summary(Sigma2RanEff)
plot(Betas)
plot(Sigma2RanEff)

########################## Visits to doctor: Longitudinal logit model ########################## 
rm(list = ls())
set.seed(12345)
Data <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/9VisitDoc.csv", sep = ",", header = TRUE, quote = "")
attach(Data)
K1 <- 7; K2 <- 2; N <- 9197
b0 <- rep(0, K1); B0 <- diag(K1)
r0 <- 5; R0 <- diag(K2)
a0 <- 0.001; d0 <- 0.001
RegLogit <- glm(DocVis ~ Age + Male + Sport + LogInc + GoodHealth + BadHealth, family = binomial(link = "logit"))
SumLogit <- summary(RegLogit)
Beta0 <- SumLogit[["coefficients"]][,1]
mcmc <- 10000; burnin <- 1000; thin <- 10
# MCMChlogit
Resultshlogit <- MCMCpack::MCMChlogit(fixed = DocVis ~ Age + Male + Sport + LogInc + GoodHealth + BadHealth, random = ~Sozh, group = "id",
                                        data = Data, burnin = burnin, mcmc = mcmc, thin = thin, 
                                        mubeta = b0, Vbeta = B0,
                                        r = r0, R = R0, nu = a0, delta = d0,
                                        beta.start = Beta0, FixOD = 0)

Betas <- Resultshlogit[["mcmc"]][,1:K1]
Sigma2RanEff <- Resultshlogit[["mcmc"]][,c(K2*N+K1+1, 2*N+K1+K2^2)]
summary(Betas)
summary(Sigma2RanEff)
plot(Betas)
plot(Sigma2RanEff)

########################## Simulation exercise: Longitudinal poisson model ########################## 
rm(list = ls())
set.seed(010101)
NT <- 1000
N <- 50
id <- c(1:N, sample(1:N, NT - N,replace=TRUE))
table(id)
x1 <- rnorm(NT); x2 <- rnorm(NT); x3 <- rnorm(NT) 
X <- cbind(1, x1, x2, x3)
K1 <- dim(X)[2]
w1 <- rnorm(NT) 
W <- cbind(1, w1)
K2 <- dim(W)[2]
B <- c(0.5, 0.4, 0.6, -0.6)
D <- c(0.7, 0.6)
sig2 <- 0.1
b1 <- rnorm(N, 0, sd = D[1]^0.5)
b2 <- rnorm(N, 0, sd = D[2]^0.5)
b <- cbind(b1, b2)
yl <- NULL
for(i in 1:NT){
  ylmeani <- X[i,]%*%B + W[i,]%*%b[id[i],]
  yli <- rnorm(1, ylmeani, sig2^0.5)
  yl <- c(yl, yli)
}
lambdait <- exp(yl)
y <- rpois(NT, lambdait)
Data <- as.data.frame(cbind(y, x1, x2, x3, w1, id))
mcmc <- 15000; burnin <- 5000; thin <- 10; tot <- mcmc + burnin
b0 <- rep(0, K1); B0 <- diag(K1); B0i <- solve(B0) 
r0 <- K2; R0 <- diag(K2); a0 <- 0.001; d0 <- 0.001
RegPois <- glm(y ~ X - 1, family = poisson(link = "log"))
SumPois <- summary(RegPois)
Beta0 <- SumPois[["coefficients"]][,1]
sig20 <- sum(SumPois[["deviance.resid"]]^2)/SumPois[["df.residual"]]
# MCMChlogit: It seems that the MCMChpoisson command is not working :-) I send an email to Professor Park reporting this error (jongheepark@snu.ac.kr) 
ResultshPoisson <- MCMCpack::MCMChpoisson(fixed = y~x1 + x2 + x3, random = ~w1, group="id",
                                          data = Data, burnin = burnin, mcmc = mcmc, thin = thin, 
                                          mubeta = b0, Vbeta = B0,
                                          r = r0, R = R0, nu = a0, delta = d0,
                                          beta.start = Beta0, FixOD = 0, sigma2.start = sig20)
Betas <- ResultshPoisson[["mcmc"]][,1:K1]
Sigma2RanEff <- ResultshPoisson[["mcmc"]][,c(K2*N+K1+1, 2*N+K1+K2^2)]
Sigma2 <- ResultshPoisson[["mcmc"]][,K2*N+K1+K2^2+1]
summary(Betas)
summary(Sigma2RanEff)
summary(Sigma2)
plot(Betas)
plot(Sigma2RanEff)
coda::geweke.diag(Betas)
coda::raftery.diag(Betas,q=0.5, r=0.05, s = 0.95)
coda::heidel.diag(Betas)

LatentMHV1 <- function(tuning, Beta, bs, sig2){
  ylhat <- rep(0, NT)
  accept <- NULL
  for(i in 1:NT){
    ids <- which(id == i)
    yi <- y[i]
    ylhatmeani <- X[i,]%*%Beta + W[i,]%*%bs[id[i],]
    ylhati <- rnorm(1, ylhatmeani, sd = sig2^0.5)
    # ylhati <- X[i,]%*%Beta + W[i,]%*%bs[id[i],] + sig2/2
    lambdahati <- exp(ylhati)
    ei <- rnorm(1, 0, sd = tuning)
    ylpropi <- ylhati + ei
    lambdapropi <- exp(ylpropi)
    logPosthati <- sum(dpois(yi, lambdahati, log = TRUE) + dnorm(ylhati, ylhatmeani, sig2^0.5, log = TRUE))
    logPostpropi <- sum(dpois(yi, lambdapropi, log = TRUE) + dnorm(ylpropi, ylhatmeani, sig2^0.5, log = TRUE))
    alphai <- min(1, exp(logPostpropi - logPosthati))
    ui <- runif(1)
    if(ui <= alphai){
      ylhati <- ylpropi
      accepti <- 1
    }else{
      ylhati <- ylhati
      accepti <- 0
    }
    ylhat[i] <- ylhati
    accept <- c(accept, accepti)
  }
  res <- list(ylhat = ylhat, accept = mean(accept))
  return(res)
}
# Beta <- B; bs <- b; tuning <- 0.1; sig2 <- sig20
# EnsLatY <- LatentMH(tuning, Beta, bs, sig2)
PostBeta <- function(D, ylhat, sig2){
  XVX <- matrix(0, K1, K1)
  XVy <- matrix(0, K1, 1)
  for(i in 1:N){
    ids <- which(id == i)
    Ti <- length(ids)
    Wi <- W[ids, ]
    Vi <- diag(Ti)*sig2 + Wi%*%D%*%t(Wi)
    ViInv <- solve(Vi)
    Xi <- X[ids, ]
    XVXi <- t(Xi)%*%ViInv%*%Xi
    XVX <- XVX + XVXi
    yi <- ylhat[ids]
    XVyi <- t(Xi)%*%ViInv%*%yi
    XVy <- XVy + XVyi
  }
  Bn <- solve(B0i + XVX)
  bn <- Bn%*%(B0i%*%b0 + XVy)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
# D <- matrix(c(0.7, 0, 0, 0.6), 2, 2); ylhat <- yl; sig2 <- sig2
# PostBeta(D, ylhat, sig2)
Postb <- function(Beta, D, ylhat, sig2){
  Di <- solve(D)
  bis <- matrix(0, N, K2)
  for(i in 1:N){
    ids <- which(id == i)
    Wi <- W[ids, ]
    Xi <- X[ids, ]
    yi <- ylhat[ids]
    Wtei <- sig2^(-1)*t(Wi)%*%(yi - Xi%*%Beta)
    Bni <- solve(sig2^(-1)*t(Wi)%*%Wi + Di)
    bni <- Bni%*%Wtei
    bi <- MASS::mvrnorm(1, bni, Bni)
    bis[i, ] <- bi
  }
  return(bis)
}
# D <- matrix(c(0.7, 0, 0, 0.6), 2, 2); ylhat <- yl; Beta <- B; sig2 <- sig2
# Postb(Beta, D, ylhat, sig2)
PostD <- function(bs){
  rn <- r0 + N
  btb <- matrix(0, K2, K2)
  for(i in 1:N){
    bsi <- bs[i, ]
    btbi <- bsi%*%t(bsi)
    btb <- btb + btbi
  }
  Rn <- d0*R0 + btb
  Sigma <- MCMCpack::riwish(v = rn, S = Rn)
  return(Sigma)
}
PostSig2 <- function(Beta, bs, ylhat){
  an <- a0 + 0.5*NT
  ete <- 0
  for(i in 1:N){
    ids <- which(id == i)
    Xi <- X[ids, ]
    yi <- ylhat[ids]
    Wi <- W[ids, ]
    ei <- yi - Xi%*%Beta - Wi%*%bs[i, ]
    etei <- t(ei)%*%ei
    ete <- ete + etei
  }
  dn <- d0 + 0.5*ete 
  sig2 <- MCMCpack::rinvgamma(1, shape = an, scale = dn)
  return(sig2)
}
PostBetas <- matrix(0, tot, K1)
PostDs <- matrix(0, tot, K2*(K2+1)/2)
Postbs <- array(0, c(N, K2, tot))
PostSig2s <- rep(0, tot)
Accepts <- rep(NULL, tot)
RegPois <- glm(y ~ X - 1, family = poisson(link = "log"))
SumPois <- summary(RegPois)
Beta <- SumPois[["coefficients"]][,1]
sig2 <- sum(SumPois[["deviance.resid"]]^2)/SumPois[["df.residual"]]
# sig2 <- 0.1
D <- diag(K2)
bs1 <- rnorm(N, 0, sd = D[1,1]^0.5)
bs2 <- rnorm(N, 0, sd = D[2,2]^0.5)
bs <- cbind(bs1, bs2)
tuning <- 0.1; ropt <- 0.44
tunepariter <- seq(round(tot/10, 0), tot, round(tot/10, 0));   l <- 1
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  LatY <- LatentMHV1(tuning = tuning, Beta = Beta, bs = bs, sig2 = sig2)
  ylhat <- LatY[["ylhat"]]
  # ylhat <- yl
  bs <- Postb(Beta = Beta, D = D, ylhat = ylhat, sig2 = sig2)
  D <- PostD(bs = bs)
  Beta <- PostBeta(D = D, ylhat = ylhat, sig2 = sig2)
  sig2 <- PostSig2(Beta = Beta, bs = bs, ylhat = ylhat)
  PostBetas[s,] <- Beta
  PostDs[s,] <- matrixcalc::vech(D)
  Postbs[, , s] <- bs
  PostSig2s[s] <- sig2
  AcceptRate <- LatY[["accept"]]
  Accepts[s] <- AcceptRate
  if(AcceptRate > ropt){
    tuning = tuning*(2-(1-AcceptRate)/(1-ropt))
  }else{
    tuning = tuning/(2-AcceptRate/ropt)
  }
  if(s == tunepariter[l]){
    # AcceptRate <- mean(Accepts[1:s])
    # if(AcceptRate > ropt){
    #   tuning = tuning*(2-(1-AcceptRate)/(1-ropt))
    # }else{
    #   tuning = tuning/(2-AcceptRate/ropt)
    # }
    print(AcceptRate)
    l <- l + 1
  }
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
Bs <- PostBetas[keep,]
Ds <- PostDs[keep,]
bs <- Postbs[, , keep]
sig2s <- PostSig2s[keep]
mean(Accepts)
summary(coda::mcmc(Bs))
plot(coda::mcmc(Bs))
summary(coda::mcmc(Ds))
plot(coda::mcmc(Ds))
summary(coda::mcmc(sig2s))
# Convergence diagnostics
coda::geweke.diag(Bs)
coda::raftery.diag(Bs, q = 0.5, r = 0.05, s = 0.95)
coda::heidel.diag(Bs)




