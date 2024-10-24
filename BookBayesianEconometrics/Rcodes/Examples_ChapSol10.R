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
    Xm <- X[, idXs]
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