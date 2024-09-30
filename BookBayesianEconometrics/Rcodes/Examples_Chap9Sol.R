########################## Effects of capital on productivity USA states: Longitudinal normal model ########################## 
rm(list = ls())
set.seed(12345)
DataGSP <- read.csv("DataApplications/8PublicCap.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(DataGSP)
K1 <- 5; K2 <- 1
N <- length(unique(id))
y <- log(gsp)
X <- cbind(1, log(pcap), log(pc), log(emp), unemp)
W <- 1
WtW <- t(W)%*%W
table(id)
T <- 17
it <- rep(1, T)
# Hyperparameters
b0 <- rep(0, K1); B0 <- diag(K1); B0i <- solve(B0)
r0 <- 5; R0 <- diag(K2)
a0 <- 0.001; d0 <- 0.001
# MCMC parameters
mcmc <- 10000
burnin <- 5000
tot <- mcmc + burnin
thin <- 1
# Gibbs functions
PostBeta <- function(sig2, D){
  Vi <- sig2*diag(T) + as.numeric(D)*it%*%t(it)
  ViInv <- solve(Vi)
  XVX <- matrix(0, K1, K1)
  XVy <- matrix(0, K1, 1)
  for(i in 1:N){
    ids <- which(id == i)
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
  Bni <- solve(sig2^(-1)*WtW + Di)
  bis <- NULL
  for(i in 1:N){
    ids <- which(id == i)
    Xi <- X[ids, ]
    yi <- y[ids]
    Wi <- it
    Wtei <- sig2^(-1)*t(Wi)%*%(yi - Xi%*%Beta)
    bni <- Bni%*%Wtei
    bi <- MASS::mvrnorm(1, bni, Bni)
    bis <- c(bis, bi)
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
    Wi <- it
    ei <- yi - Xi%*%Beta - bs[i]*Wi
    etei <- t(ei)%*%ei
    ete <- ete + etei
  }
  dn <- d0 + 0.5*ete 
  sig2 <- invgamma::rinvgamma(1, shape = an, rate = dn)
  return(sig2)
}
PostD <- function(bs){
  rn <- r0 + N
  btb <- 0
  for(i in 1:N){
    bsi <- bs[i]
    btbi <- bsi%*%t(bsi)
    btb <- btb + btbi
  }
  Rn <- r0*R0 + btb
  Sigma <- LaplacesDemon::rinvwishart(nu = rn, S = Rn)
  return(Sigma)
}
PostBetas <- matrix(0, tot, K1)
PostDs <- matrix(0, tot, K2*(K2+1)/2)
PostSig2s <- rep(0, tot)
Postbs <- matrix(0, tot, N)
RegLS <- lm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp)
SumLS <- summary(RegLS)
Beta <- SumLS[["coefficients"]][,1]
sig2 <- SumLS[["sigma"]]^2
D <- diag(K2)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  bs <- Postb(Beta = Beta, sig2 = sig2, D = D)
  D <- PostD(bs = bs)
  # D <- matrix(0.1,1,1)
  Beta <- PostBeta(sig2 = sig2, D = D)
  sig2 <- PostSig2(Beta = Beta, bs = bs, D = D)
  # sig2 <- 1.453e-03
  PostBetas[s,] <- Beta
  PostDs[s,] <- matrixcalc::vech(D)
  PostSig2s[s] <- sig2
  Postbs[s,] <- bs
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
Bs <- PostBetas[keep,]
Ds <- PostDs[keep,]
bs <- Postbs[keep,]
sig2s <- PostSig2s[keep]
summary(coda::mcmc(Bs))
summary(coda::mcmc(Ds))
summary(coda::mcmc(bs))
summary(coda::mcmc(sig2s))
summary(coda::mcmc(Ds^0.5/(sig2s+Ds)^0.5))
# Convergence diagnostics
coda::geweke.diag(Bs)
coda::raftery.diag(Bs,q=0.5,r=0.05,s = 0.95)
coda::heidel.diag(Bs)
