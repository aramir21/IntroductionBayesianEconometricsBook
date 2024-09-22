########################## Effects of institutions: Multivariate regression ########################## 
rm(list = ls())
set.seed(010101)
DataInst <- read.csv("DataApplications/4Institutions.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(DataInst)
Y <- cbind(logpcGDP95, PAER)
X <- cbind(1, logMort, Africa, Asia, Other)
M <- dim(Y)[2]
K <- dim(X)[2]
# Hyperparameters
B0 <- matrix(0, K, M)
c0 <- 100
V0 <- c0*diag(K)
Psi0 <- 5*diag(M)
a0 <- 5
S <- 10000 #Number of posterior draws
betadraw = matrix(double(S*K*M), ncol=K*M)
Sigmadraw = matrix(double(S*M*M), ncol=M*M)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
for (s in 1:S) {
  Results <- bayesm::rmultireg(Y, X, Bbar = B0, A = solve(V0), nu = a0, V = Psi0)
  betadraw[s,] <- Results$B
  Sigmadraw[s,] <- Results$Sigma
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0), "% done"))
}
close(pb)
summary(coda::mcmc(betadraw))
summary(coda::mcmc(Sigmadraw))
########################## Demand and supply: Simulation ##########################
# Simulation
rm(list = ls())
set.seed(12345)
B1 <- 5; B2 <- -0.5; B3 <- 0.8; B4 <- -0.4; B5 <- 0.7; SD <- 0.5
A1 <- -2; A2 <- 0.5; A3 <- -0.4; SS <- 0.5
P0 <- (A1-B1)/(B2-A2); P2 <- -B3/(B2-A2); P3 <- -B4/(B2-A2); P1 <- A3/(B2-A2); P4 <- -B5/(B2-A2)
T0 <- B1+B2*P0; T2 <- B3+B2*P2; T3 <- B4+B2*P3; T1 <- B2*P1; T4 <- B5+B2*P4;
n <- 5000
ED <- rnorm(n, 0, SD); ES <- rnorm(n, 0, SS)
VP <- (ES-ED)/(B2-A2); UQ <- B2*VP+ED
y <- rnorm(n, 10, 1); pc <- rnorm(n, 5, 1); er <- rnorm(n, 15, 1); ps <- rnorm(n, 5, 1);
p <- P0+P1*er+P2*y+P3*pc+P4*ps+VP
q <- T0+T1*er+T2*y+T3*pc+T4*ps+UQ
#Inference
Y <- cbind(p, q)
X <- cbind(1, er, y, pc, ps)
M <- dim(Y)[2]
K <- dim(X)[2]
# Hyperparameters
B0 <- matrix(0, K, M)
c0 <- 100
V0 <- c0*diag(K)
Psi0 <- 5*diag(M)
a0 <- 5
S <- 10000 #Number of posterior draws
betadraw = matrix(double(S*K*M), ncol=K*M)
Sigmadraw = matrix(double(S*M*M), ncol=M*M)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
for (s in 1:S) {
  Results <- bayesm::rmultireg(Y, X, Bbar = B0, A = solve(V0), nu = a0, V = Psi0)
  betadraw[s,] <- Results$B
  Sigmadraw[s,] <- Results$Sigma
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0), "% done"))
}
close(pb)
summary(coda::mcmc(betadraw))
summary(coda::mcmc(Sigmadraw))
beta2 <- betadraw[,7]/betadraw[,2] # Effect of price on demand 
summary(coda::mcmc(beta2))
beta3 <- betadraw[,8] - beta2*betadraw[,3] # Effect of income on demand
summary(coda::mcmc(beta3))
beta4 <- betadraw[,9] - beta2*betadraw[,4] # Effect of complementary price on demand
summary(coda::mcmc(beta4))
beta5 <- betadraw[,10] - beta2*betadraw[,5] # Effect of substitute price  on demand
summary(coda::mcmc(beta5))
# Parameters of the supply are overidentified, that is, there are different ways to get them
alpha3_1 <- -beta3*betadraw[,2]/betadraw[,3] # Effect of substitute price  on demand
summary(coda::mcmc(alpha3_1))
alpha3_2 <- -beta4*betadraw[,2]/betadraw[,4] # Effect of substitute price  on demand
summary(coda::mcmc(alpha3_2))
alpha2_2 <- beta2 + beta3/betadraw[,3] # Effect of substitute price  on demand
summary(coda::mcmc(alpha2_2))
alpha2_3 <- beta2 + beta4/betadraw[,4] # Effect of substitute price  on demand
summary(coda::mcmc(alpha2_3))

########################## SUR: Utilities demand ########################## 
rm(list = ls())
set.seed(010101)
library(dplyr)
DataUt <- read.csv("DataApplications/Utilities.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
DataUtEst <- DataUt %>%  
  filter(Electricity != 0 & Water !=0 & Gas != 0)
attach(DataUtEst)
y1 <- log(Electricity); y2 <- log(Water); y3 <- log(Gas)
X1 <- cbind(1, LnPriceElect, LnPriceWater, LnPriceGas, IndSocio1, IndSocio2, Altitude, Nrooms, HouseholdMem, Lnincome)
X2 <- cbind(1, LnPriceElect, LnPriceWater, LnPriceGas, IndSocio1, IndSocio2, Nrooms, HouseholdMem)
X3 <- cbind(1, LnPriceElect, LnPriceWater, LnPriceGas, IndSocio1, IndSocio2, Altitude, Nrooms, HouseholdMem)
y <- c(y1, y2, y3)
X <- as.matrix(Matrix::bdiag(X1, X2, X3))
M <- 3; K1 <- dim(X1)[2]; K2 <- dim(X2)[2]; K3 <- dim(X3)[2] 
K <- K1 + K2 + K3
N <- length(y1)
# Hyperparameters
b0 <- rep(0, K)
c0 <- 100
B0 <- c0*diag(K)
B0i <- solve(B0)
Psi0 <- 5*diag(M)
Psi0i <- solve(Psi0)
a0 <- M
IN <- diag(N)
an <- a0 + N
#Posterior draws
S <- 6000 #Number of posterior draws
burnin <- 1000
thin <- 1
tot <- S+burnin
# Gibbs functions
PostBeta <- function(Sigma){
  Aux <- solve(Sigma)%x%IN
  Bn <- solve(B0i + t(X)%*%Aux%*%X)
  bn <- Bn%*%(B0i%*%b0 + t(X)%*%Aux%*%y)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostSigma <- function(Beta){
  B1 <- Beta[1:K1]; B2 <- Beta[(K1+1):(K1+K2)]; B3 <- Beta[(K1+K2+1):(K1+K2+K3)]
  U1 <- y1 - X1%*%B1; U2 <- y2 - X2%*%B2; U3 <- y3 - X3%*%B3
  U <- cbind(U1, U2, U3)
  Psin <- solve(Psi0i + t(U)%*%U)
  Sigmai <- rWishart::rWishart(1, df = an, Sigma = Psin)
  Sigma <- solve(Sigmai[,,1]) 
  return(Sigma)
}
PostBetas <- matrix(0, tot, K)
PostSigmas <- matrix(0, tot, M*(M+1)/2)
Beta <- rep(1, K)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Sigma <- PostSigma(Beta = Beta)
  Beta <- PostBeta(Sigma = Sigma)
  PostBetas[s,] <- Beta
  PostSigmas[s,] <- matrixcalc::vech(Sigma)
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
Bs <- PostBetas[keep,]
Names <- c("Const", "LnPriceElect", "LnPriceWater", "LnPriceGas", "IndSocio1", "IndSocio2", 
           "Altitude", "Nrooms", "HouseholdMem", "Lnincome", "Const",
           "LnPriceElect", "LnPriceWater", "LnPriceGas", "IndSocio1", "IndSocio2", 
           "Nrooms", "HouseholdMem","Const",
           "LnPriceElect", "LnPriceWater", "LnPriceGas", "IndSocio1", "IndSocio2", 
           "Altitude", "Nrooms", "HouseholdMem")
colnames(Bs) <- Names
summary(coda::mcmc(Bs))
Sigmas <- PostSigmas[keep,]
summary(coda::mcmc(Sigmas))
# Gibbs sampling 2 using sums rather than Kronecker product
PostBetaNew <- function(Sigma){
  Sigmai <- solve(Sigma)
  XSX <- matrix(0, K, K)
  XSy <- matrix(0, K, 1)
  for(i in 1:N){
    Xi <- as.matrix(Matrix::bdiag(t(X1[i,]), t(X2[i,]), t(X3[i,])))
    yi <- c(y1[i], y2[i], y3[i])
    XSXi <- t(Xi)%*%Sigmai%*%Xi
    XSyi <- t(Xi)%*%Sigmai%*%yi
    XSX <- XSX + XSXi
    XSy <- XSy + XSyi
  }
  Bn <- solve(B0i + XSX)
  bn <- Bn%*%(B0i%*%b0 + XSy)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostBetas <- matrix(0, tot, K)
PostSigmas <- matrix(0, tot, M*(M+1)/2)
Beta <- rep(1, K)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Sigma <- PostSigma(Beta = Beta)
  Beta <- PostBetaNew(Sigma = Sigma)
  PostBetas[s,] <- Beta
  PostSigmas[s,] <- matrixcalc::vech(Sigma)
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
Bs <- PostBetas[keep,]
Names <- c("Const", "LnPriceElect", "LnPriceWater", "LnPriceGas", "IndSocio1", "IndSocio2", 
           "Altitude", "Nrooms", "HouseholdMem", "Lnincome", "Const",
           "LnPriceElect", "LnPriceWater", "LnPriceGas", "IndSocio1", "IndSocio2", 
           "Nrooms", "HouseholdMem","Const",
           "LnPriceElect", "LnPriceWater", "LnPriceGas", "IndSocio1", "IndSocio2", 
           "Altitude", "Nrooms", "HouseholdMem")
colnames(Bs) <- Names
summary(coda::mcmc(Bs))
Sigmas <- PostSigmas[keep,]
summary(coda::mcmc(Sigmas))

########################## Instrumental variables: Simulation ########################## 
# Weak instruments
rm(list = ls())
set.seed(010101)
N <- 100
k <- 2
B <- rep(1, k)
G <- c(1, 0.2)
s12 <- 0.8
SIGMA <- matrix(c(1, s12, s12, 1), 2, 2)
z <- rnorm(N); # w <- rnorm(N)
Z <- cbind(1, z); w <- matrix(1,N,1)
S <- 100
U <- replicate(S, MASS::mvrnorm(n = N, mu = rep(0, 2), SIGMA))
x <- G[1] + G[2]*z + U[,2,]
y <- B[1] + B[2]*x + U[,1,]
VarX <- G[2]^2+1 # Population variance of x
EU1U2 <- s12 # Covariance U1
BiasPopB2 <- EU1U2/VarX
# Hyperparameters
d0 <- 0.001/2
a0 <- 0.001/2
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
g0 <- rep(0, 2)
G0 <- 1000*diag(2)
G0i <- solve(G0)
nu <- 3
Psi0 <- nu*diag(2)
# MCMC parameters
mcmc <- 5000
burnin <- 1000
tot <- mcmc + burnin
thin <- 1
# Gibbs sampling
Gibbs <- function(x, y){
  Data <- list(y = y, x = x, w = w, z = Z)
  Mcmc <- list(R = mcmc, keep = thin, nprint = 0)
  Prior <- list(md = g0, Ad = G0i, mbg = b0, Abg = B0i, nu = nu, V = Psi0)
  RestIV <- bayesm::rivGibbs(Data = Data, Mcmc = Mcmc, Prior = Prior)
  PostBIV <- mean(RestIV[["betadraw"]])
  ResLM <- MCMCpack::MCMCregress(y ~ x + w - 1, b0 = b0, B0 = B0i, c0 = a0, d0 = d0)
  PostB <- mean(ResLM[,1])
  Res <- c(PostB,PostBIV)
  return(Res)
}
PosteriorMeans <- sapply(1:S, function(s) {Gibbs(x = x[,s], y = y[,s])})
rowMeans(PosteriorMeans)
Model <- c(replicate(S, "Ordinary"), replicate(S, "Instrumental"))
postmeans <- c(t(PosteriorMeans))
df <- data.frame(postmeans, Model, stringsAsFactors = FALSE)
library(ggplot2); library(latex2exp)
histExo <- ggplot(df, aes(x = postmeans, fill = Model)) +
  geom_histogram(bins = 40, position = "identity", color = "black", alpha = 0.5) +
  labs(title = "Overlayed Histograms", x = "Value", y = "Count") +
  scale_fill_manual(values = c("blue", "red")) +
  geom_vline(aes(xintercept = mean(postmeans[1:S])), color = "black", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(postmeans[101:200])), color = "black", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = B[2]), color = "green", linewidth = 1, linetype = "dashed") +
  xlab(TeX("$E[\\beta_2]$")) + ylab("Frequency") + 
  ggtitle("Histogram: Posterior means simulating 100 samples") 
histExo  

# Endogenous instruments
rm(list = ls())
set.seed(010101)
N <- 100
k <- 2
B <- rep(1, k)
G <- c(1, 1)
delta <- 0.2 # Effect of instrument on y
s12 <- 0.8
SIGMA <- matrix(c(1, s12, s12, 1), 2, 2)
z <- rnorm(N); # w <- rnorm(N)
Z <- cbind(1, z); w <- matrix(1,N,1)
S <- 100
U <- replicate(S, MASS::mvrnorm(n = N, mu = rep(0, 2), SIGMA))
x <- G[1] + G[2]*z + U[,2,]
y <- B[1] + B[2]*x + delta*z + U[,1,]
VarX <- G[2]^2+1 # Population variance of x
EU1U2 <- s12 # Covariance U1
BiasPopB2 <- EU1U2/VarX
# Hyperparameters
d0 <- 0.001/2
a0 <- 0.001/2
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
g0 <- rep(0, 2)
G0 <- 1000*diag(2)
G0i <- solve(G0)
nu <- 3
Psi0 <- nu*diag(2)
# MCMC parameters
mcmc <- 5000
burnin <- 1000
tot <- mcmc + burnin
thin <- 1
# Gibbs sampling
Gibbs <- function(x, y){
  Data <- list(y = y, x = x, w = w, z = Z)
  Mcmc <- list(R = mcmc, keep = thin, nprint = 0)
  Prior <- list(md = g0, Ad = G0i, mbg = b0, Abg = B0i, nu = nu, V = Psi0)
  RestIV <- bayesm::rivGibbs(Data = Data, Mcmc = Mcmc, Prior = Prior)
  PostBIV <- mean(RestIV[["betadraw"]])
  ResLM <- MCMCpack::MCMCregress(y ~ x + w - 1, b0 = b0, B0 = B0i, c0 = a0, d0 = d0)
  PostB <- mean(ResLM[,1])
  Res <- c(PostB,PostBIV)
  return(Res)
}
PosteriorMeans <- sapply(1:S, function(s) {Gibbs(x = x[,s], y = y[,s])})
rowMeans(PosteriorMeans)
Model <- c(replicate(S, "Ordinary"), replicate(S, "Instrumental"))
postmeans <- c(t(PosteriorMeans))
df <- data.frame(postmeans, Model, stringsAsFactors = FALSE)
library(ggplot2); library(latex2exp)
histExo <- ggplot(df, aes(x = postmeans, fill = Model)) +
  geom_histogram(bins = 40, position = "identity", color = "black", alpha = 0.5) +
  labs(title = "Overlayed Histograms", x = "Value", y = "Count") +
  scale_fill_manual(values = c("blue", "red")) +
  geom_vline(aes(xintercept = mean(postmeans[1:S])), color = "black", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = mean(postmeans[101:200])), color = "black", linewidth = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = B[2]), color = "green", linewidth = 1, linetype = "dashed") +
  xlab(TeX("$E[\\beta_2]$")) + ylab("Frequency") + 
  ggtitle("Histogram: Posterior means simulating 100 samples") 
histExo  

########################## Instrumental variables: Simulation and Gibbs sampler ########################## 
rm(list = ls())
set.seed(010101)
N <- 100
k <- 2
kz <- 2
B <- rep(1, k)
G <- rep(1, kz)
s12 <- 0.8
SIGMA <- matrix(c(1, s12, s12, 1), 2, 2)
z <- rnorm(N);
Z <- cbind(1, z); w <- matrix(1,N,1)
S <- 100
U <- MASS::mvrnorm(n = N, mu = rep(0, 2), SIGMA)
x <- G[1] + G[2]*z + U[,2]
y <- B[1] + B[2]*x + U[,1]
X <- cbind(w, x)
# Hyperparameters
d0 <- 0.001/2
a0 <- 0.001/2
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
g0 <- rep(0, kz)
G0 <- 1000*diag(kz)
G0i <- solve(G0)
nu <- 3
Psi0 <- nu*diag(2)
Psi0i <- solve(Psi0)
# MCMC parameters
mcmc <- 5000
burnin <- 1000
tot <- mcmc + burnin
thin <- 1
# Auxiliary elements
XtX <- t(X)%*%X 
ZtZ <- t(Z)%*%Z 
nun <- nu + N
# Gibbs sampling
PostBeta <- function(Sigma, Gamma){
  w1 <- Sigma[1,1] - Sigma[1,2]^2/Sigma[2,2]
  Bn <- solve(w1^(-1)*XtX + B0i)
  yaux <- y - (Sigma[1,2]/Sigma[2,2])*(x - Z%*%Gamma)
  bn <- Bn%*%(B0i%*%b0 + w1^(-1)*t(X)%*%yaux)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostGamma <- function(Sigma, Beta){
  w2 <- Sigma[2,2] - Sigma[1,2]^2/Sigma[1,1]
  Gn <- solve(w2^(-1)*ZtZ + G0i)
  xaux <- x - (Sigma[1,2]/Sigma[1,1])*(y - X%*%Beta)
  gn <- Gn%*%(G0i%*%g0 + w2^(-1)*t(Z)%*%xaux)
  Gamma <- MASS::mvrnorm(1, gn, Gn)
  return(Gamma)
}
PostSigma <- function(Beta, Gamma){
  Uy <- y - X%*%Beta; Ux <- x - Z%*%Gamma
  U <- cbind(Uy, Ux)
  Psin <- solve(Psi0i + t(U)%*%U)
  Sigmai <- rWishart::rWishart(1, df = nun, Sigma = Psin)
  Sigma <- solve(Sigmai[,,1]) 
  return(Sigma)
}
PostBetas <- matrix(0, tot, k)
PostGammas <- matrix(0, tot, kz)
PostSigmas <- matrix(0, tot, 2*(2+1)/2)
Beta <- rep(0, k); Gamma <- rep(0, kz)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Sigma <- PostSigma(Beta = Beta, Gamma = Gamma)
  Beta <- PostBeta(Sigma = Sigma, Gamma = Gamma)
  Gamma <- PostGamma(Sigma = Sigma, Beta = Beta)
  PostBetas[s,] <- Beta
  PostGammas[s,] <- Gamma
  PostSigmas[s,] <- matrixcalc::vech(Sigma)
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
Bs <- PostBetas[keep,]
Gs <- PostGammas[keep,]
Sigmas <- PostSigmas[keep,]
summary(coda::mcmc(Bs))
summary(coda::mcmc(Gs))
summary(coda::mcmc(Sigmas))

########################## Instrumental variables: Effect of institutions on GDP ########################## 
# Weak instruments
rm(list = ls())
set.seed(010101)
DataInst <- read.csv("DataApplications/6Institutions.csv", sep = ",", header = TRUE, fileEncoding = "latin1")
attach(DataInst)
y <- logpcGDP95; x <- PAER
w <- cbind(1, Africa, Asia, Other); Z <- cbind(1, logMort)
# Hyperparameters
k <- 5; kz <- 2
b0 <- rep(0, k)
c0 <- 100
B0 <- c0*diag(k)
B0i <- solve(B0)
g0 <- rep(0, kz)
G0 <- 100*diag(kz)
G0i <- solve(G0)
nu <- 5
Psi0 <- nu*diag(2)
# MCMC parameters
mcmc <- 50000
burnin <- 10000
tot <- mcmc + burnin
thin <- 5
# Gibbs sampling
Data <- list(y = y, x = x, w = w, z = Z)
Mcmc <- list(R = mcmc, keep = thin, nprint = 100)
Prior <- list(md = g0, Ad = G0i, mbg = b0, Abg = B0i, nu = nu, V = Psi0)
RestIV <- bayesm::rivGibbs(Data = Data, Mcmc = Mcmc, Prior = Prior)
summary(RestIV[["deltadraw"]])
summary(coda::mcmc(RestIV[["betadraw"]]))
summary(RestIV[["gammadraw"]])
summary(RestIV[["Sigmadraw"]])

########################## Multivariate probit: Simulation ########################## 
remove(list = ls())
n<-500  #Individuals
p<-2  #Number dependent variables (For instance to buy or not to buy different products)
xo<- 2  #Number of choice dependent variables (For instance price of products)
xi<- 1  #Number of regressors that dependt on individuals (For instance income)
B1<-c(0.5,-1.2,0.7,0.8)
B2<-c(1.5,-0.8,0.5,0) # The last regressor is not relevant for the second product
Cor<-matrix(c(1,0.5,0.5,1),p,p)
yX<-NULL
for (i in 1:n){
  ei<-MASS::mvrnorm(1,c(0,0),Cor)
  xs<-rnorm(xi)
  x1i<-c(1,rnorm(xo),xs)
  yl1i<-sum(x1i*B1)
  y1i<-yl1i+ei[1]>0
  x2i<-c(1,rnorm(xo),xs)
  yl2i<-sum(x2i*B2)
  y2i<-yl2i+ei[2]>0
  yXi<-rbind(c(y1i,x1i),c(y2i,x2i))
  yXi<-cbind(i,yXi)
  colnames(yXi)<-c("id","y","cte","x1o","x2o","xi")
  yX<-rbind(yX,yXi)
}

nd <- 4 # Number of regressors
Xd <- as.matrix(yX[,3:6])
XcreateMP<-function(p,nxs,nind,Data){
  pandterm = function(message) {
    stop(message, call. = FALSE)
  }
  if (missing(nxs)) 
    pandterm("requires number of regressors: include intercept if required")
  if (missing(nind)) 
    pandterm("requires number of units (individuals)")
  if (missing(Data)) 
    pandterm("requires dataset")
  if (nrow(Data)!=nind*2)
    pandterm("check dataset! number of units times number alternatives should be equal to dataset rows")
  
  XXDat<-array(0,c(p,1+nxs,nind))
  XX<-array(0,c(p,nxs*p,nind))
  YY<-array(0,c(p,1,nind))
  is<- seq(p,nind*p,p)
  cis<- seq(nxs,nxs*p+1,nxs)
  for(i in is){
    j<-which(i==is)
    XXDat[,,j]<-as.matrix(Data[c((i-(p-1)):i),-1])
    YY[,,j]<-XXDat[,1,j]
    for(l in 1:p){
      XX[l,((cis[l]-(nxs-1)):cis[l]),j]<-XXDat[l,-1,j]
    }
  }
  return(list(y=YY,X=XX))
}
Dat <- XcreateMP(p = p, nxs = nd, nind = n, Data = yX)
y<-NULL
X<-NULL
for(i in 1:dim(Dat$y)[3]){
  y<-c(y,Dat$y[,,i])
  X<-rbind(X,Dat$X[,,i])
}
DataMP = list(p=p, y=y, X=X)
# Hyperparameters
k <- dim(X)[2]
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)
B0i <- solve(B0)
a0 <- p - 1 + 3
Psi0 <- a0*diag(p)
Prior <- list(betabar = b0, A = B0i, nu = a0, V = Psi0)
# MCMC parameters
mcmc <- 2000
thin <- 5
Mcmc <- list(R = mcmc, keep = thin)
Results <- bayesm::rmvpGibbs(Data = DataMP, Mcmc = Mcmc, Prior = Prior)
betatildeEq1 <- Results$betadraw[,1:4] / sqrt(Results$sigmadraw[,1])
summary(coda::mcmc(betatildeEq1))
betatildeEq2 <- Results$betadraw[,5:8] / sqrt(Results$sigmadraw[,4])
summary(coda::mcmc(betatildeEq2))
sigmadraw12 <-  Results$sigmadraw[,2] / (Results$sigmadraw[,1]*Results$sigmadraw[,4])^0.5
summary(coda::mcmc(sigmadraw12))

### Gibs sampler using basically the SUR structure augmenting with latent variables yl1 and yl2
ids1 <- seq(1, p*n, 2); ids2 <- seq(2, p*n, 2)  
y1 <- y[ids1]; y2 <- y[ids2]
X1 <- X[ids1,1:4]; X2 <- X[ids2, 5:7]  
y <- c(y1, y2)
X <- as.matrix(Matrix::bdiag(X1, X2))
M <- 2; K1 <- dim(X1)[2]; K2 <- dim(X2)[2]
K <- K1 + K2
N <- length(y1)
# Hyperparameters
b0 <- rep(0, K)
c0 <- 1000
B0 <- c0*diag(K)
B0i <- solve(B0)
Psi0 <- 4*diag(M)
Psi0i <- solve(Psi0)
a0 <- M
IN <- diag(N)
an <- a0 + N
#Posterior draws
S <- 2000 #Number of posterior draws
burnin <- 1
thin <- 5
tot <- S+burnin
# Gibbs functions
PostBeta <- function(Sigma, Yl1, Yl2){
  Aux <- solve(Sigma)%x%IN
  Bn <- solve(B0i + t(X)%*%Aux%*%X)
  Yl <- c(Yl1, Yl2)
  bn <- Bn%*%(B0i%*%b0 + t(X)%*%Aux%*%Yl)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostSigma <- function(Beta, Yl1, Yl2){
  B1 <- Beta[1:K1]; B2 <- Beta[(K1+1):(K1+K2)]
  U1 <- Yl1 - X1%*%B1; U2 <- Yl2 - X2%*%B2
  U <- cbind(U1, U2)
  Psin <- solve(Psi0i + t(U)%*%U)
  Sigmai <- rWishart::rWishart(1, df = an, Sigma = Psin)
  Sigma <- solve(Sigmai[,,1]) 
  return(Sigma)
}
PostYl1 <- function(Beta, Sigma, i){
  Beta1 <- Beta[1:K1]
  Ylmean1 <- X1[i,]%*%Beta1
  if(y1[i] == 1){
    Yl1i <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = Ylmean1, sd = 1)
  }else{
    Yl1i <- truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = Ylmean1, sd = 1)
  }
  return(Yl1i)
}
PostYl2 <- function(Beta, Sigma, i){
  Beta2 <- Beta[(K1+1):(K1+K2)]
  Ylmean2 <- X2[i,]%*%Beta2
  if(y2[i] == 1){
    Yl2i <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = Ylmean2, sd = 1)
  }else{
    Yl2i <- truncnorm::rtruncnorm(1, a = -Inf, b = 0, mean = Ylmean2, sd = 1)
  }
  return(Yl2i)
}

PostBetas <- matrix(0, tot, K)
PostSigmas <- matrix(0, tot, M*(M+1)/2)
Beta <- rep(1, K)
Sigma <- diag(M)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Yl1 <- sapply(1:n, function(i){PostYl1(Beta = Beta, Sigma = Sigma, i)})
  Yl2 <- sapply(1:n, function(i){PostYl2(Beta = Beta, Sigma = Sigma, i)})
  Sigma <- PostSigma(Beta = Beta, Yl1 = Yl1, Yl2 = Yl2)
  Beta <- PostBeta(Sigma = Sigma, Yl1 = Yl1, Yl2 = Yl2)
  PostBetas[s,] <- Beta
  PostSigmas[s,] <- matrixcalc::vech(Sigma)
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
Bs <- PostBetas[keep,]
betatilde1 <- Bs[,1:K1] / sqrt(PostSigmas[keep,1])
summary(coda::mcmc(betatilde1))
betatilde2 <- Bs[,(K1+1):(K1+K2)] / sqrt(PostSigmas[keep,3])
summary(coda::mcmc(betatilde2))
sigmadraw12 <-  PostSigmas[keep,2] / (PostSigmas[keep,1]*PostSigmas[keep,3])^0.5
summary(coda::mcmc(sigmadraw12))

### Gibs sampler using setting in Section 8.4 in book





