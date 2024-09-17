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




