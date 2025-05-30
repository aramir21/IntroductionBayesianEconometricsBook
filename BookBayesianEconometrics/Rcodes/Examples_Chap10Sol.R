########################## Simulation exercise: BMA normal model unique models ########################## 
rm(list = ls())
set.seed(010101)
N <- 1000
K1 <- 6; K2 <- 4; K3 <- 30; K <- K1 + K2 + K3
X1 <- matrix(rnorm(N*K1,1 ,1), N, K1)
X2 <- matrix(rbinom(N*K2, 1, 0.5), N, K2)
X3 <- matrix(rnorm(N*K3,1 ,1), N, K3)
X <- cbind(X1, X2, X3)
e <- rnorm(N, 0, 0.5)
B <- c(1,0,0,0,0.5,0,0,0,0,-0.7, rep(0, 30))
y <- 1 + X%*%B + e
#### Markov chain Monte Carlo model composition from scratch
LogMLfunt <- function(Model){
  indr <- Model == 1
  kr <- sum(indr)
  if(kr > 0){
    gr <- ifelse(N > kr^2, 1/N, kr^(-2))
    Xr <- matrix(Xnew[ , indr], ncol = kr)
    # PX <- diag(N) - Xr%*%solve(t(Xr)%*%Xr)%*%t(Xr)
    # s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
    PX <- Xr%*%solve(t(Xr)%*%Xr)%*%t(Xr)
    s2pos <- c((t(y - mean(y))%*%(y - mean(y))) - t(y)%*%PX%*%y/(1 + gr))
    mllMod <- (kr/2)*log(gr/(1+gr))-(N-1)/2*log(s2pos)
  }else{
    gr <- ifelse(N > kr^2, 1/N, kr^(-2))
    # PX <- diag(N)
    # s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
    s2pos <- c((t(y - mean(y))%*%(y - mean(y))))
    mllMod <- (kr/2)*log(gr/(1+gr))-(N-1)/2*log(s2pos)
  }
  return(mllMod)
}
Xnew <- apply(X, 2, scale)
# MC3 method
# Initial models
M <- 100
Models <- matrix(rbinom(K*M, 1, p = 0.5), ncol = K, nrow = M + 800)
Models <- unique(Models)[1:M,]
mllnew <- sapply(1:M, function(s){LogMLfunt(matrix(Models[s,], 1, K))})
oind <- order(mllnew, decreasing = TRUE)
mllnew <- mllnew[oind]
Models <- Models[oind, ]
# Hyperparameters MC3
iter <- 50000
pb <- winProgressBar(title = "progress bar", min = 0, max = iter, width = 300)
s <- 1
while(s <= iter){
  ActModel <- Models[M,]
  idK <- which(ActModel == 1)
  Kact <- length(idK)
  Continue <- 0
  while(Continue == 0){
    if(Kact < K & Kact > 1){
      CardMol <- K
      opt <- sample(1:3, 1)
      if(opt == 1){ # Same
        CandModel <- ActModel
      }else{
        if(opt == 2){ # Add
          All <- 1:K
          NewX <- sample(All[-idK], 1)
          CandModel <- ActModel
          CandModel[NewX] <- 1
        }else{ # Subtract
          LessX <- sample(idK, 1)
          CandModel <- ActModel
          CandModel[LessX] <- 0
        }
      }
    }else{
      CardMol <- K + 1
      if(Kact == K){
        opt <- sample(1:2, 1)
        if(opt == 1){ # Same
          CandModel <- ActModel
        }else{ # Subtract
          LessX <- sample(1:K, 1)
          CandModel <- ActModel
          CandModel[LessX] <- 0
        }
      }else{
        if(K == 1){
          opt <- sample(1:3, 1)
          if(opt == 1){ # Same
            CandModel <- ActModel
          }else{
            if(opt == 2){ # Add
              All <- 1:K
              NewX <- sample(All[-idK], 1)
              CandModel <- ActModel
              CandModel[NewX] <- 1
            }else{ # Subtract
              LessX <- sample(idK, 1)
              CandModel <- ActModel
              CandModel[LessX] <- 0
            }
          }
        }else{ # Add
          NewX <- sample(1:K, 1)
          CandModel <- ActModel
          CandModel[NewX] <- 1
        }
      }
    }
    check <- NULL
    for(j in 1:M){
      if(sum(Models[j,] == CandModel) == K){
        checkj <- 0
        check <- c(check, checkj)
      }else{
          checkj <- 1
          check <- c(check, checkj)
        }
      }
    dimUniModels <- sum(check)
    if(dimUniModels == M){
      Continue <- 1
    }else{
      Continue <- 0
    }
  }
  LogMLact <- LogMLfunt(matrix(ActModel, 1, K))
  LogMLcand <- LogMLfunt(matrix(CandModel, 1, K))
  alpha <- min(1, exp(LogMLcand-LogMLact)) # Let's reasonably assume same prior model probability for candidate and actual, and same carnality of neighbor models
  u <- runif(1)
  if(u <= alpha){
    mllnew[M] <- LogMLcand
    Models[M, ] <- CandModel
    oind <- order(mllnew, decreasing = TRUE)
    mllnew <- mllnew[oind]
    Models <- Models[oind, ]
  }else{
    mllnew <- mllnew
    Models <- Models
  }
  s <- s + 1
  setWinProgressBar(pb, s, title=paste( round(s/iter*100, 0),"% done"))
}
close(pb)
ModelsUni <- unique(Models)
mllnewUni <- sapply(1:dim(ModelsUni)[1], function(s){LogMLfunt(matrix(ModelsUni[s,], 1, K))})
StMarLik <- exp(mllnewUni-mllnewUni[1])
PMP <- StMarLik/sum(StMarLik) # PMP based on unique selected models
plot(PMP)
PIP <- NULL
for(k in 1:K){
  PIPk <- sum(PMP[which(ModelsUni[,k] == 1)])
  PIP <- c(PIP, PIPk)
}
plot(PIP)
Means <- matrix(0, M, K)
Vars <- matrix(0, M, K)
for(m in 1:M){
  idXs <- which(ModelsUni[m,] == 1)
  if(length(idXs) == 0){
    Regm <- lm(y ~ 1)
  }else{
    Xm <- X[, idXs]
    Regm <- lm(y ~ Xm)
    SumRegm <- summary(Regm)
    Means[m, idXs] <- SumRegm[["coefficients"]][-1,1]
    Vars[m, idXs] <- SumRegm[["coefficients"]][-1,2]^2 
  }
}
BMAmeans <- colSums(Means*PMP)
BMAsd <- (colSums(PMP*Vars)  + colSums(PMP*(Means-matrix(rep(BMAmeans, each = M), M, K))^2))^0.5 
plot(BMAmeans)
plot(BMAsd)
plot(BMAmeans/BMAsd)
########################## Determinants of export diversification: BMA normal model ########################## 
rm(list = ls())
set.seed(010101)
Data <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/10ExportDiversificationHHI.csv", sep = ",", header = TRUE, quote = "")
attach(Data)
y <- Data[,1]
X <- as.matrix(Data[,-1])
K <- dim(X)[2]
#### BIC approximation
# BMAglm <- BMA::bic.glm(avghhi ~ island+landlocked+opec+british+dutch+french+spanish+portuguese+language+commonlaw+europe+asia+africa+namerica+       
#                        smamerica+avgpolity+avggov+avggoveff+avgtradefreedom+avglinfl+avgedu1+avgedu2+avgedu3+avgedu4+avgedu5+avgexports+avgimports+
#                        avgnatres+avgcapital+avgfdi+avglgdpcap+avgtot+avgfuel+avgki+avglpop+avgroads, data = Data, glm.family = gaussian(), strict = FALSE, OR = 50) 
BMAglm <- BMA::bicreg(X, y, strict = FALSE, OR = 50)
summary(BMAglm)

#### Markov chain Monte Carlo model composition using BMA package
BMAreg <- BMA::MC3.REG(y, X, num.its=10000)
Models <- unique(BMAreg[["variables"]])
nModels <- dim(Models)[1]
nVistModels <- dim(BMAreg[["variables"]])[1]
PMPmc3 <- NULL
for(m in 1:nModels){
  idModm <- NULL
  for(j in 1:nVistModels){
    if(sum(Models[m,] == BMAreg[["variables"]][j,]) == K){
      idModm <- c(idModm, j)
    }else{
      idModm <- idModm
    } 
  }
  PMPm <- sum(BMAreg[["post.prob"]][idModm])
  PMPmc3 <- c(PMPmc3, PMPm)
}
PMPmc3
PIPmc3 <- NULL
for(k in 1:K){
  PIPk <- sum(PMPmc3[which(Models[,k] == 1)])
  PIPmc3 <- c(PIPmc3, PIPk)
}
plot(PIPmc3)
Meansmc3 <- matrix(0, nModels, K)
Varsmc3 <- matrix(0, nModels, K)
for(m in 1:nModels){
  idXs <- which(Models[m,] == 1)
  if(length(idXs) == 0){
    Regm <- lm(y ~ 1)
  }else{
    Xm <- X[, idXs]
    Regm <- lm(y ~ Xm)
    SumRegm <- summary(Regm)
    Meansmc3[m, idXs] <- SumRegm[["coefficients"]][-1,1]
    Varsmc3[m, idXs] <- SumRegm[["coefficients"]][-1,2]^2 
  }
}
BMAmeansmc3 <- colSums(Meansmc3*PMPmc3)
BMAsdmc3 <- (colSums(PMPmc3*Varsmc3)  + colSums(PMPmc3*(Meansmc3-matrix(rep(BMAmeansmc3, each = nModels), nModels, K))^2))^0.5 
plot(BMAmeansmc3)
plot(BMAsdmc3)
plot(BMAmeansmc3/BMAsdmc3)
Ratio <- BMAmeansmc3/BMAsdmc3
RessulMC3 <- as.data.frame(cbind(PIPmc3, Models[1,], BMAmeansmc3, BMAsdmc3, Ratio))
PMPmc3
########################## Simulation exercise: IV BMA ########################## 
rm(list = ls())
set.seed(010101)
simIV <- function(delta1,delta2,beta0,betas1,betas2,beta2,Sigma,n,z) {
  eps <- matrix(rnorm(3*n),ncol=3) %*% chol(Sigma)
  xs1 <- z%*%delta1 + eps[,1]
  xs2 <- z%*%delta2 + eps[,2]
  x2 <- rnorm(dim(z)[1])
  y <- beta0+betas1*xs1+betas2*xs2+beta2*x2 + eps[,3]
  X <- as.matrix(cbind(xs1,xs2,1,x2)) 
  colnames(X) <- c("x1en","x2en","cte","xex")
  y <- matrix(y,dim(z)[1],1)
  colnames(y) <- c("y")
  list(X=X,y=y)
}
n <- 1000 ; p <- 3 
z <- matrix(runif(n*p),ncol=p)
rho31 <- 0.8; rho32 <- 0.5;
Sigma <- matrix(c(1,0,rho31,0,1,rho32,rho31,rho32,1),ncol=3)
delta1 <- c(4,-1,2); delta2 <- c(-2,3,-1); betas1 <- .5; betas2 <- -1; beta2 <- 1; beta0 <- 2
simiv <- simIV(delta1,delta2,beta0,betas1,betas2,beta2,Sigma,n,z)
nW <- 18
W <- matrix(rnorm(nW*dim(z)[1]),dim(z)[1],nW)
YXW<-cbind(simiv$y, simiv$X, W)
y <- YXW[,1]; X <- YXW[,2:3]; W <- YXW[,-c(1:4)]
Xnew <- cbind(X, W)
BMAglm <- BMA::bicreg(Xnew, y, strict = FALSE, OR = 50) 
summary(BMAglm)
BMAreg <- BMA::MC3.REG(y, Xnew, num.its=10000)
Models <- unique(BMAreg[["variables"]])
nModels <- dim(Models)[1]
nVistModels <- dim(BMAreg[["variables"]])[1]
K <- dim(Xnew)[2]
PMP <- NULL
for(m in 1:nModels){
  idModm <- NULL
  for(j in 1:nVistModels){
    if(sum(Models[m,] == BMAreg[["variables"]][j,]) == K){
      idModm <- c(idModm, j)
    }else{
      idModm <- idModm
    } 
  }
  PMPm <- sum(BMAreg[["post.prob"]][idModm])
  PMP <- c(PMP, PMPm)
}
PMP
PIP <- NULL
for(k in 1:K){
  PIPk <- sum(PMP[which(Models[,k] == 1)])
  PIP <- c(PIP, PIPk)
}
plot(PIP)
Means <- matrix(0, nModels, K)
Vars <- matrix(0, nModels, K)
for(m in 1:nModels){
  idXs <- which(Models[m,] == 1)
  if(length(idXs) == 0){
    Regm <- lm(y ~ 1)
  }else{
    Xm <- Xnew[, idXs]
    Regm <- lm(y ~ Xm)
    SumRegm <- summary(Regm)
    Means[m, idXs] <- SumRegm[["coefficients"]][-1,1]
    Vars[m, idXs] <- SumRegm[["coefficients"]][-1,2]^2 
  }
}
BMAmeans <- colSums(Means*PMP)
BMAsd <- (colSums(PMP*Vars)  + colSums(PMP*(Means-matrix(rep(BMAmeans, each = nModels), nModels, K))^2))^0.5 
plot(BMAmeans)
plot(BMAsd)
plot(BMAmeans/BMAsd)

#PIP from BIC, and coefficients from Instrumental variable
BMAglm <- BMA::bicreg(Xnew, y, strict = FALSE, OR = 50) 
summary(BMAglm)
PMPBIC <- BMAglm[["postprob"]]
ModelsBIC <- BMAglm[["which"]]
nModels <- dim(ModelsBIC)[1]
K <- dim(Xnew)[2]
Means <- matrix(0, nModels, K)
Vars <- matrix(0, nModels, K)
for(m in 1:nModels){
  idXs <- which(ModelsBIC[m,] == 1)
  if(length(idXs) == 0){
    Regm <- lm(y ~ 1)
  }else{
    Xm <- Xnew[, idXs]
    Regm <- ivreg::ivreg(y ~ Xm | z + W)
    SumRegm <- summary(Regm)
    Means[m, idXs] <- SumRegm[["coefficients"]][-1,1]
    Vars[m, idXs] <- SumRegm[["coefficients"]][-1,2]^2 
  }
}
BMAmeans <- colSums(Means*PMPBIC)
BMAsd <- (colSums(PMPBIC*Vars)  + colSums(PMPBIC*(Means-matrix(rep(BMAmeans, each = nModels), nModels, K))^2))^0.5 
plot(BMAmeans)
plot(BMAsd)
plot(BMAmeans/BMAsd)

########################## Determinants of export diversification: BMA normal model ########################## 
rm(list = ls())
set.seed(010101)
DataMain <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/11ExportDiversificationHHI.csv", sep = ",", header = TRUE, quote = "")
DataInst <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/12ExportDiversificationHHIInstr.csv", sep = ",", header = TRUE, quote = "")
attach(DataMain)
attach(DataInst)
y <- DataMain[,1]
X <- as.matrix(DataMain[,2])
W <- as.matrix(DataMain[,-c(1:2)])
Z <- as.matrix(DataInst)
S <- 10000; burnin <- 1000
regivBMA <- ivbma::ivbma(Y = y, X = X, Z = Z, W = W, s = S+burnin, b = burnin, odens = S, print.every = round(S/10), run.diagnostics = FALSE)
PIPmain <- regivBMA[["L.bar"]] # PIP outcome
PIPmain
EVmain <- regivBMA[["rho.bar"]] # Posterior mean outcome
EVmain
PIPaux <- regivBMA[["M.bar"]] # PIP auxiliary
PIPaux
EVaux <- regivBMA[["lambda.bar"]] # Posterior mean auxiliary
plot(EVaux[,1])
EVsigma <- regivBMA[["Sigma.bar"]] # Posterior mean variance matrix
EVsigma
summary(coda::mcmc(regivBMA[["Sigma"]][1,2,]))

########################## Determinants of internet adoption: BMA logit model ########################## 
rm(list = ls())
set.seed(010101)
DataMain <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/13InternetMed.csv", sep = ",", header = TRUE, quote = "")
attach(DataMain)
### BIC approximation
y <- DataMain[,1]
X <- as.matrix(DataMain[,-1])
BMAglm <- BMA::bic.glm(X, y, strict = FALSE, OR = 50, glm.family = binomial(link="logit"))
summary(BMAglm)

########################## Determinants of market value of soccer players in Europe: BMA gamma model ########################## 
rm(list = ls())
set.seed(010101)
DataMain <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/14ValueFootballPlayers.csv", sep = ",", header = TRUE, quote = "")
attach(DataMain)
### BIC approximation
y <- DataMain[,1]
X <- as.matrix(DataMain[,-1])
BMAglm <- BMA::bic.glm(X, y, strict = FALSE, OR = 50, glm.family = Gamma(link="log"))
summary(BMAglm)

########################## Determinants of fertility in Bostwana: BMA Poisson model ########################## 
rm(list = ls())
set.seed(010101)
DataMain <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/15Fertil2.csv", sep = ",", header = TRUE, quote = "")
attach(DataMain)
### BIC approximation
y <- DataMain[,1]
X <- as.matrix(DataMain[,-1])
BMAglm <- BMA::bic.glm(X, y, strict = FALSE, OR = 50, glm.family = poisson(link="log"))
summary(BMAglm)

########################## Simulation exercise: BMA Logit ########################## 
rm(list = ls())
set.seed(010101)
n<-1000 
B<-c(0.5,0.8,-1.2)
X<-matrix(cbind(rep(1,n),rnorm(n,0,1),rnorm(n,0,1)),n,length(B))
p <- exp(X%*%B)/(1+exp(X%*%B))
y <- rbinom(n, 1, p)
table(y)
nXgar<-25
Xgar<-matrix(rnorm(nXgar*n),n,nXgar)
df<-as.data.frame(cbind(y,X[,-1],Xgar))
colnames(df) <- c("y", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27")
### BIC approximation
BMAglmLogit <- BMA::bic.glm(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27, data = df, glm.family = binomial(link="logit"), strict = FALSE, OR = 50)
summary(BMAglmLogit)

########################## Simulation exercise: BMA logit using BIC from scratch ########################## 
Xnew <- df[,-1]
BICfunt <- function(Model){
  indr <- Model == 1
  kr <- sum(indr)
  if(kr > 0){
    Xr <- as.matrix(Xnew[ , indr])
    model <- glm(y ~ Xr, family = binomial(link = "logit"))
    model_bic <- BIC(model)
    mllMod <- -model_bic/2
  }else{
    model <- glm(y ~ 1, family = binomial(link = "logit"))
    model_bic <- BIC(model)
    mllMod <- -model_bic/2
  }
  return(mllMod)
}
# MC3 method
# Initial models
M <- 500
K <- dim(df)[2] - 1
Models <- matrix(rbinom(K*M, 1, p = 0.5), ncol = K, nrow = M)
mllnew <- sapply(1:M, function(s){BICfunt(matrix(Models[s,], 1, K))})
oind <- order(mllnew, decreasing = TRUE)
mllnew <- mllnew[oind]
Models <- Models[oind, ]
# Hyperparameters MC3
iter <- 25000
pb <- winProgressBar(title = "progress bar", min = 0, max = iter, width = 300)
s <- 1
while(s <= iter){
  ActModel <- Models[M,]
  idK <- which(ActModel == 1)
  Kact <- length(idK)
  if(Kact < K & Kact > 1){
    CardMol <- K
    opt <- sample(1:3, 1)
    if(opt == 1){ # Same
      CandModel <- ActModel
    }else{
      if(opt == 2){ # Add
        All <- 1:K
        NewX <- sample(All[-idK], 1)
        CandModel <- ActModel
        CandModel[NewX] <- 1
      }else{ # Subtract
        LessX <- sample(idK, 1)
        CandModel <- ActModel
        CandModel[LessX] <- 0
      }
    }
  }else{
    CardMol <- K + 1
    if(Kact == K){
      opt <- sample(1:2, 1)
      if(opt == 1){ # Same
        CandModel <- ActModel
      }else{ # Subtract
        LessX <- sample(1:K, 1)
        CandModel <- ActModel
        CandModel[LessX] <- 0
      }
    }else{
      if(K == 1){
        opt <- sample(1:3, 1)
        if(opt == 1){ # Same
          CandModel <- ActModel
        }else{
          if(opt == 2){ # Add
            All <- 1:K
            NewX <- sample(All[-idK], 1)
            CandModel <- ActModel
            CandModel[NewX] <- 1
          }else{ # Subtract
            LessX <- sample(idK, 1)
            CandModel <- ActModel
            CandModel[LessX] <- 0
          }
        }
      }else{ # Add
        NewX <- sample(1:K, 1)
        CandModel <- ActModel
        CandModel[NewX] <- 1
      }
    }
  }
  LogMLact <- BICfunt(matrix(ActModel, 1, K))
  LogMLcand <- BICfunt(matrix(CandModel, 1, K))
  alpha <- min(1, exp(LogMLcand-LogMLact)) # Let's reasonably assume same prior model probability for candidate and actual, and same carnality of neighbor models
  u <- runif(1)
  if(u <= alpha){
    mllnew[M] <- LogMLcand
    Models[M, ] <- CandModel
    oind <- order(mllnew, decreasing = TRUE)
    mllnew <- mllnew[oind]
    Models <- Models[oind, ]
  }else{
    mllnew <- mllnew
    Models <- Models
  }
  s <- s + 1
  setWinProgressBar(pb, s, title=paste( round(s/iter*100, 0),"% done"))
}
close(pb)
ModelsUni <- unique(Models)
mllnewUni <- sapply(1:dim(ModelsUni)[1], function(s){BICfunt(matrix(ModelsUni[s,], 1, K))})
StMarLik <- exp(mllnewUni-mllnewUni[1])
PMP <- StMarLik/sum(StMarLik) # PMP based on unique selected models
plot(PMP)
ModelsUni[1,]
PIP <- NULL
for(k in 1:K){
  PIPk <- sum(PMP[which(ModelsUni[,k] == 1)])
  PIP <- c(PIP, PIPk)
}
plot(PIP)
nModels <- dim(ModelsUni)[1]
Means <- matrix(0, nModels, K)
Vars <- matrix(0, nModels, K)
for(m in 1:nModels){
  idXs <- which(ModelsUni[m,] == 1)
  if(length(idXs) == 0){
    Regm <- glm(y ~ 1, family = binomial(link = "logit"))
  }else{
    Xm <- as.matrix(Xnew[, idXs])
    Regm <- glm(y ~ Xm, family = binomial(link = "logit"))
    SumRegm <- summary(Regm)
    Means[m, idXs] <- SumRegm[["coefficients"]][-1,1]
    Vars[m, idXs] <- SumRegm[["coefficients"]][-1,2]^2 
  }
}
BMAmeans <- colSums(Means*PMP)
BMAsd <- (colSums(PMP*Vars)  + colSums(PMP*(Means-matrix(rep(BMAmeans, each = nModels), nModels, K))^2))^0.5 
plot(BMAmeans)
plot(BMAsd)
plot(BMAmeans/BMAsd)

########################## Dynamic Bayesian model average: Application Exchange rate COP USD ########################## 
rm(list = ls())
set.seed(010101)
DataMain <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/19ExchangeRateCOPUSD.csv", sep = ",", header = TRUE, quote = "")
X <- DataMain[-c(1:12,227),c(4, 5, 7, 9, 11, 13, 14, 15)]
DifIntRat <- X[,8] - X[,7]
DifInf <- X[,1] - X[,2]
DifGrowth <- X[,3] - X[,4]
DifM3 <- X[,6] - X[,5]
XDif <- cbind(DifIntRat, DifInf, DifGrowth, DifM3)
matplot(XDif, type = "l")
y <- DataMain[-c(1:13),17] # Var exchange rate depends on lag of other variables
plot(y, type = "l")
IRP <- c(1, 0, 0, 0); PPP <- c(0, 1, 0, 0)
Taylor <- c(0, 1, 1, 0); Money <- c(0, 0, 1, 1)
Models <- rbind(IRP, PPP, Taylor, Money)
T0 <- 50
dma.test <- dma::dma(XDif, y, Models, lambda=.99, gamma=.99, initialperiod = T0)
plot(ts(dma.test[["pmp"]][,1], start = 2006, end = 2023, frequency = 12), type = "l", main = "Exchange rate models", ylab = "Posterior model probability", ylim = c(0,1))
lines(ts(dma.test[["pmp"]][,2], start = 2006, end = 2023, frequency = 12), col = "red")
lines(ts(dma.test[["pmp"]][,3], start = 2006, end = 2023, frequency = 12), col = "green")
lines(ts(dma.test[["pmp"]][,4], start = 2006, end = 2023, frequency = 12), col = "blue")
legend(x = 2009, y = 0.6, legend = c("IRP", "PPP", "Taylor rule", "Money supply"), col = c("black", "red", "green", "blue"), lty=1:1, cex=0.8)
require(latex2exp)
plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5,
                fill = "blue") +
    geom_line(aes(y = mean), colour = "black", linewidth = 0.8) +
    geom_hline(yintercept=0, colour = "red") +
    ylab(TeX("$\\beta_{Growth_t}^{Money}$")) + xlab("Time") 
  print(p)
}
menathetaIRP <- dma.test[["thetahat"]][4,,4]
SDthetaIRP <- (dma.test[["Vtheta"]][4,,4])^0.5
LimInf <- menathetaIRP - 2*SDthetaIRP
LimSup <- menathetaIRP + 2*SDthetaIRP
# Time <- DataMain$date[-c(1:12,215)] 
library(ggplot2,tibble,dplyr)
df <- tibble(t = 1:214,
             mean = menathetaIRP,
             lower = LimInf,
             upper = LimSup)
Fig <- plot_filtering_estimates(df)

########################## Dynamic Bayesian model average in logistic model: Simulation exercise ########################## 
rm(list = ls())
set.seed(010101)
T <- 1100; K <- 3
X <- matrix(rnorm(T*K, mean = 0.5, sd = 0.8), T, K)
combs <- expand.grid(c(0,1), c(0,1), c(0,1))
B1 <- 0.5
B2t <- seq(1, 2, length.out=T )
a <- 0.5
B3t <- c(rep(-1,round(a*T)), rep(0,round((1-a)*T)))
B4 <- 1.2
yl <- B1 + X[,1]*B2t + X[,2]*B3t + X[,3]*B4
py <- exp(yl)/(1 + exp(yl))
y <- rbinom(T, 1, prob = py)
table(y)
T0 <- 100
dma.test <- dma::logistic.dma(X, y, combs[-1,], lambda = 0.99, alpha = 0.99, initialsamp = T0)
plot(dma.test[["pmp"]][7, -c(1:T0)], type = "l", col = "green", main = "Posterior model probability: Model all regressors vs model regressors 1 and 3", xlab = "Time", ylab = "PMP", ylim = c(0, 1))
lines(dma.test[["pmp"]][5,-c(1:T0)], col = "red")
legend(x = 0, y = 0.9, legend = c("Model: All regressors", "Model: Regressors 1 and 3"), col = c("green", "red"), lty=1:1, cex=0.8)
dma.test1 <- dma::logistic.dma(X, y, combs[-1,], lambda = 0.99, alpha = 0.95, initialsamp = T0)
plot(dma.test1[["pmp"]][7, -c(1:T0)], type = "l", col = "green", main = "Posterior model probability: Model all regressors vs model regressors 1 and 3", xlab = "Time", ylab = "PMP", ylim = c(0, 1))
lines(dma.test1[["pmp"]][5,-c(1:T0)], col = "red")
legend(x = 0, y = 0.9, legend = c("Model: All regressors", "Model: Regressors 1 and 3"), col = c("green", "red"), lty=1:1, cex=0.8)
require(latex2exp)
plot(dma.test[["theta"]][7,-c(1:T0),1], type = "l", col = "green", main = "Bayesian model average filtering recursion", xlab = "Time", ylab = TeX("$\\beta_{1}$"))
abline(h = B1, col = "red")
legend(x = 0, y = 0.4, legend = c("State filtering", "State population"), col = c("green", "red"), lty=1:1, cex=0.8)
plot(dma.test[["theta"]][7,-c(1:T0),2], type = "l", col = "green", main = "Bayesian model average filtering recursion", xlab = "Time", ylab = TeX("$\\beta_{2t}$"), ylim = c(0.5,2))
lines(B2t[-c(1:T0)], col = "red")
legend(x = 0, y = 0.8, legend = c("State filtering", "State population"), col = c("green", "red"), lty=1:1, cex=0.8)
plot(dma.test[["theta"]][7,-c(1:T0),3], type = "l", col = "green", main = "Bayesian model average filtering recursion", xlab = "Time", ylab = TeX("$\\beta_{3t}$"), ylim = c(-1.4,0))
lines(B3t[-c(1:T0)], col = "red")
legend(x = 0, y = -0.4, legend = c("State filtering", "State population"), col = c("green", "red"), lty=1:1, cex=0.8)
plot(dma.test[["theta"]][7,-c(1:T0),4], type = "l", col = "green", main = "Bayesian model average filtering recursion", xlab = "Time", ylab = TeX("$\\beta_{4t}$"))
abline(h = B4, col = "red")
legend(x = 0, y = 1.3, legend = c("State filtering", "State population"), col = c("green", "red"), lty=1:1, cex=0.8)

