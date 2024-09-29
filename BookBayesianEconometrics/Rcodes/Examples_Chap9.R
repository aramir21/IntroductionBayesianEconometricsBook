########################## Effects of capital on productivity USA states: Longitudinal normal model ########################## 
rm(list = ls())
set.seed(12345)
DataGSP <- read.csv("DataApplications/8PublicCap.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
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
set.seed(12345)
NT <- 1000
N <- 25
id <- c(1:N, sample(1:N, NT - N,replace=TRUE))
table(id)
x1 <- rnorm(NT); x2 <- rnorm(NT); x3 <- rnorm(NT) 
X <- cbind(1, x1, x2, x3)
K1 <- dim(X)[2]
w1 <- rnorm(NT) 
W <- cbind(1, w1)
K2 <- dim(W)[2]
B <- c(0.5, 0.4, 0.6, -0.6)
D <- c(0.5, 0.6)
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
r0 <- 5; R0 <- diag(K2); a0 <- 0.001; d0 <- 0.001
# MCMChregress
Resultshreg <- MCMCpack::MCMChregress(fixed = y~x1 + x2 + x3, random = ~w1, group="id",
                      data = Data, burnin = burnin, mcmc = mcmc, thin = thin, 
                      mubeta = b0, Vbeta = B0,
                      r = r0, R = R0, nu = a0, delta = d0)
Betas <- Resultshreg[["mcmc"]][,1:K1]
Sigma2RanEff <- Resultshreg[["mcmc"]][,c(55, 58)]
Sigma2 <- Resultshreg[["mcmc"]][,59]
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
PostSig2 <- function(Beta, D, bs){
  an <- a0 + 0.5*N*T
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
  dn <- 1/d0 + 0.5*ete 
  sig2 <- invgamma::rinvgamma(1, shape = an, scale = dn)
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
  Rn <- R0 + btb
  Sigma <- LaplacesDemon::rinvwishart(nu = rn, S = Rn)
  return(Sigma)
}
PostBetas <- matrix(0, tot, K1)
PostDs <- matrix(0, tot, K2*(K2+1)/2)
PostSig2s <- rep(0, tot)
Postbs <- array(0, c(N, K2, tot))
Beta <- rep(1, K1)
D <- diag(K2)
RegLS <- lm(y ~ X - 1)
SumLS <- summary(RegLS)
sig2 <- SumLS[["sigma"]]^0.5
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  bs <- Postb(Beta = Beta, sig2 = sig2, D = D)
  D <- PostD(bs = bs)
  Beta <- PostBeta(sig2 = sig2, D = D)
  sig2 <- PostSig2(Beta = Beta, bs = bs, D = D)
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
summary(coda::mcmc(Ds))
summary(coda::mcmc(sig2s))
# Convergence diagnostics
coda::geweke.diag(Bs)
coda::raftery.diag(Bs,q=0.5,r=0.05,s = 0.95)
coda::heidel.diag(Bs)


