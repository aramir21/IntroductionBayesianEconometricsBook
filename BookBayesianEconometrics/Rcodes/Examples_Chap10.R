########################## Simulation exercise: BMA normal model ########################## 
rm(list = ls())
set.seed(010101)
N <- 1000
K1 <- 6; K2 <- 4; K <- K1 + K2
X1 <- matrix(rnorm(N*K1,1 ,1), N, K1)
X2 <- matrix(rbinom(N*K2, 1, 0.5), N, K2)
X <- cbind(X1, X2)
e <- rnorm(N, 0, 0.5)
B <- c(1,0,0,0,0.5,0,0,0,0,-0.7)
y <- 1 + X%*%B + e
#### BIC approximation
BMAglm <- BMA::bicreg(X, y, strict = FALSE, OR = 50) 
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
combs <- expand.grid(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1),c(0,1), c(0,1), c(0,1), c(0,1), c(0,1))
Xnew <- apply(X, 2, scale)
mll <- sapply(1:2^K, function(s){LogMLfunt(matrix(combs[s,], 1, K))})
# Results
MaxPMP <- which.max(mll)
StMarLik <- exp(mll-max(mll))
PMP <- StMarLik/sum(StMarLik)
PMP[MaxPMP]
combs[MaxPMP,]
PIP <- NULL
for(k in 1:K){
  PIPk <- sum(PMP[which(combs[,k] == 1)])
  PIP <- c(PIP, PIPk)
}
plot(PIP)
nModels <- dim(combs)[1]
Means <- matrix(0, nModels, K)
Vars <- matrix(0, nModels, K)
for(m in 1:nModels){
  idXs <- which(combs[m,] == 1)
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
# MC3 method
# Initial models
M <- 100
Models <- matrix(rbinom(K*M, 1, p = 0.5), ncol = K, nrow = M)
mllnew <- sapply(1:M, function(s){LogMLfunt(matrix(Models[s,], 1, K))})
oind <- order(mllnew, decreasing = TRUE)
mllnew <- mllnew[oind]
Models <- Models[oind, ]
# Hyperparameters MC3
iter <- 6000
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
# PMP using number of visits
nModels <- dim(ModelsUni)[1]
StMarLik <- exp(mllnew-mllnew[1])
PMPold <- StMarLik/sum(StMarLik) # PMP all selected models
PMPot <- NULL
PMPap <- NULL
FreqMod <- NULL
for(m in 1:nModels){
  idModm <- NULL
  for(j in 1:M){
    if(sum(ModelsUni[m,] == Models[j,]) == K){
      idModm <- c(idModm, j)
    }else{
      idModm <- idModm
    }
  }
  PMPm <- sum(PMPold[idModm]) # PMP unique models using sum of all selected models
  PMPot <- c(PMPot, PMPm)
  PMPapm <- length(idModm)/M # PMP using relative frequency in all selected models
  PMPap <- c(PMPap, PMPapm)
  FreqMod <- c(FreqMod, length(idModm))
}
cbind(PMP, PMPot, PMPap)
PIP <- NULL
for(k in 1:K){
  PIPk <- sum(PMP[which(ModelsUni[,k] == 1)])
  PIP <- c(PIP, PIPk)
}
plot(PIP)
Means <- matrix(0, nModels, K)
Vars <- matrix(0, nModels, K)
for(m in 1:nModels){
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
BMAsd <- (colSums(PMP*Vars)  + colSums(PMP*(Means-matrix(rep(BMAmeans, each = nModels), nModels, K))^2))^0.5 
plot(BMAmeans)
plot(BMAsd)
plot(BMAmeans/BMAsd)

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
y <- YXW[,1]; X <- YXW[,2:3]; W <- YXW[,-c(1:3)]
S <- 10000; burnin <- 1000
regivBMA <- ivbma::ivbma(Y = y, X = X, Z = z, W = W, s = S+burnin, b = burnin, odens = S, print.every = round(S/10), run.diagnostics = FALSE)
PIPmain <- regivBMA[["L.bar"]] # PIP outcome
EVmain <- regivBMA[["rho.bar"]] # Posterior mean outcome
plot(EVmain)
PIPaux <- regivBMA[["M.bar"]] # PIP auxiliary
EVaux <- regivBMA[["lambda.bar"]] # Posterior mean auxiliary
plot(EVaux[,1])
plot(EVaux[,2])
EVsigma <- regivBMA[["Sigma.bar"]] # Posterior mean variance matrix
EVsigma 

########################## Simulation exercise: BMA logit ########################## 
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

########################## Simulation exercise: BMA gamma ########################## 
rm(list = ls())
set.seed(010101)
n<-1000 
B<- c(0.5, 0.2, 0.1) # c(0.5,0.5,0.5)
X<-matrix(cbind(rep(1,n),rnorm(n,0,0.5),rnorm(n,0,0.5)),n,length(B))
y1 <- (X%*%B)^(-1) # exp(X%*%B)
plot(y1)
y <- rgamma(n,y1,scale=1)
nXgar<-25
Xgar<-matrix(rnorm(nXgar*n),n,nXgar)
df<-as.data.frame(cbind(y,X[,-1],Xgar))
colnames(df) <- c("y", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27")
### BIC approximation
BMAglmGamma <- BMA::bic.glm(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27, data = df, glm.family = Gamma(link="inverse"), strict = FALSE, OR = 50)
summary(BMAglmGamma)

########################## Simulation exercise: BMA Poisson ########################## 
rm(list = ls())
set.seed(010101)
n<-1000 
B<-c(2,1.1,0.7)
X<-matrix(cbind(rep(1,n),rnorm(n,0,1),rnorm(n,0,1)),n,length(B))
y1<-exp(X%*%B)
y<-rpois(n,y1)
nXgar<-25
Xgar<-matrix(rnorm(nXgar*n),n,nXgar)
df<-as.data.frame(cbind(y,X[,-1],Xgar))
colnames(df) <- c("y", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27")
### BIC approximation
BMAglmPoisson <- BMA::bic.glm(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21+x22+x23+x24+x25+x26+x27, data = df, glm.family = poisson(link="log"), strict = FALSE, OR = 50)
summary(BMAglmPoisson)

########################## Simulation exercise: BMA Poisson using BIC from scratch ########################## 
Xnew <- df[,-1]
BICfunt <- function(Model){
  indr <- Model == 1
  kr <- sum(indr)
  if(kr > 0){
    Xr <- as.matrix(Xnew[ , indr])
    model <- glm(y ~ Xr, family = poisson(link = "log"))
    model_bic <- BIC(model)
    mllMod <- -model_bic/2
  }else{
    model <- glm(y ~ 1, family = poisson(link = "log"))
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
    Regm <- glm(y ~ 1, family = poisson(link = "log"))
  }else{
    Xm <- as.matrix(Xnew[, idXs])
    Regm <- glm(y ~ Xm, family = poisson(link = "log"))
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





