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
iter <- 10000
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

df <- as.data.frame(cbind(y, X)) 
colnames(df) <- c("y", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30", "x31", "x32", "x33", "x34", "x35", "x36", "x37", "x38", "x39", "x40")
write.csv(df, file="51SimNormalBMANew.csv", row.names=FALSE)
#### BIC approximation
BMAglm <- BMA::bic.glm(y ~ X, data = df, glm.family = gaussian()) 
summary(BMAglm)

#### Markov chain Monte Carlo model composition using BMA package
BMAreg <- BMA::MC3.REG(y, X, num.its=500)
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




