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
Xnew <- apply(df[,-1], 2, scale)
# Xnew <- df[,-1]
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
Xnew <- df[,-1]
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

########################## Dynamic Bayesian model average: Simulation exercise ########################## 
rm(list = ls())
set.seed(010101)
T <- 500; K <- 3
X <- matrix(rnorm(T*K, mean = 0.5, sd = 0.8), T, K)
combs <- expand.grid(c(0,1), c(0,1), c(0,1))
B1 <- 0.5
B2t <- seq(1, 2, length.out=T )
a <- 0.75
B3t <- c(rep(-1,round(a*T)), rep(0,round((1-a)*T)))
B4 <- 1.2
sigma <- 1; mu <- rnorm(T, 0, sigma)
y <- B1 + X[,1]*B2t + X[,2]*B3t + X[,3]*B4 + mu
T0 <- 50
dma.test <- dma::dma(X, y, combs, lambda=.99, gamma=.99, initialperiod = T0)
plot(dma.test[["pmp"]][-c(1:T0),8], type = "l", col = "green", main = "Posterior model probability: Model all regressors vs model regressors 2 and 4", xlab = "Time", ylab = "PMP")
lines(dma.test[["pmp"]][-c(1:T0),6], col = "red")
legend(x = 0, y = 0.9, legend = c("Model: All regressors", "Model: Regressors 2 and 4"), col = c("green", "red"), lty=1:1, cex=0.8)
require(latex2exp)
plot(dma.test[["thetahat.ma"]][-c(1:T0),1], type = "l", col = "green", main = "Bayesian model average filtering recursion", xlab = "Time", ylab = TeX("$\\beta_{1}$"))
abline(h = B1, col = "red")
legend(x = 0, y = 0.4, legend = c("State filtering", "State population"), col = c("green", "red"), lty=1:1, cex=0.8)
plot(dma.test[["thetahat.ma"]][-c(1:T0),2], type = "l", col = "green", main = "Bayesian model average filtering recursion", xlab = "Time", ylab = TeX("$\\beta_{2t}$"), ylim = c(0.5,2))
lines(B2t[-c(1:T0)], col = "red")
legend(x = 0, y = 0.8, legend = c("State filtering", "State population"), col = c("green", "red"), lty=1:1, cex=0.8)
plot(dma.test[["thetahat.ma"]][-c(1:T0),3], type = "l", col = "green", main = "Bayesian model average filtering recursion", xlab = "Time", ylab = TeX("$\\beta_{3t}$"))
lines(B3t[-c(1:T0)], col = "red")
legend(x = 0, y = -0.4, legend = c("State filtering", "State population"), col = c("green", "red"), lty=1:1, cex=0.8)
plot(dma.test[["thetahat.ma"]][-c(1:T0),4], type = "l", col = "green", main = "Bayesian model average filtering recursion", xlab = "Time", ylab = TeX("$\\beta_{4t}$"))
abline(h = B4, col = "red")
legend(x = 0, y = 1.3, legend = c("State filtering", "State population"), col = c("green", "red"), lty=1:1, cex=0.8)

########################## Simulation exercise: Bayes factors ########################## 
rm(list = ls())
set.seed(010101)
N <- 500; K <- 5; K2 <- 3 
B <- c(0.7, 0.3, 0.7, -0.2, 0.2) 
X1 <- rbinom(N, 1, 0.3)
X2 <- matrix(rnorm(K2*N), N, K2)
X <- cbind(1, X1, X2)
Y <- X%*%B + rnorm(N, 0, sd = 1)
# Hyperparameters
d0 <- 4
a0 <- 4
b0 <- rep(0, K)
cOpt <- 0.5
an <- N + a0; B0 <- cOpt*diag(K)
Bn <- solve(solve(B0)+t(X)%*%X); bhat <- solve(t(X)%*%X)%*%t(X)%*%Y
bn <- Bn%*%(solve(B0)%*%b0+t(X)%*%X%*%bhat)
dn <- as.numeric(d0 + t(Y-X%*%bhat)%*%(Y-X%*%bhat)+t(bhat - b0)%*%solve(solve(t(X)%*%X)+B0)%*%(bhat - b0))
Hn <- as.matrix(Matrix::forceSymmetric(dn*Bn/an))
S <- 10000
LogMarLikLM <- function(X, c0){
  K <- dim(X)[2]
  N <- dim(X)[1]	
  # Hyperparameters
  B0 <- c0*diag(K)
  b0 <- rep(0, K)
  # Posterior parameters
  bhat <- solve(t(X)%*%X)%*%t(X)%*%Y
  # Force this matrix to be symmetric
  Bn <- as.matrix(Matrix::forceSymmetric(solve(solve(B0) + t(X)%*%X))) 
  bn <- Bn%*%(solve(B0)%*%b0 + t(X)%*%X%*%bhat)
  dn <- as.numeric(d0 + t(Y)%*%Y+t(b0)%*%solve(B0)%*%b0-t(bn)%*%solve(Bn)%*%bn)
  an <- a0 + N
  # Log marginal likelihood
  logpy <- (N/2)*log(1/pi)+(a0/2)*log(d0)-(an/2)*log(dn) + 0.5*log(det(Bn)/det(B0)) + lgamma(an/2)-lgamma(a0/2)
  return(-logpy)
}
# Empirical Bayes: Obtain c0 maximizing the log 
# EB2 <- optim(c0, fn = LogMarLikLM, method = "Brent", lower = 0.0001, upper = 10^6, X = X)
# cOpt <- EB2$par # c0 <- EB2$par
LogMarM2 <- -LogMarLikLM(X = X, c0 = cOpt)
LogMarM1 <- -LogMarLikLM(X = X[,1:4], c0 = cOpt)
BF12 <- exp(LogMarM1-LogMarM2) 
BF12; 1/BF12
2*log(1/BF12)
# Savage-Dickey density ratio
# Posterior evaluation
Brest <- 0
an <- N + a0
B0 <- cOpt*diag(K) 
Bn <- solve(solve(B0)+t(X)%*%X)
bhat <- solve(t(X)%*%X)%*%t(X)%*%Y
bn <- Bn%*%(solve(B0)%*%b0+t(X)%*%X%*%bhat)
dn <- as.numeric(d0 + t(Y-X%*%bhat)%*%(Y-X%*%bhat)+t(bhat - b0)%*%solve(solve(t(X)%*%X)+B0)%*%(bhat - b0))
Hn <- as.matrix(Matrix::forceSymmetric(dn*Bn/an))
# VarBeta5 <- Hn[5,5] - Hn[5,1:4]%*%solve(Hn[1:4,1:4])%*%Hn[1:4,5] 
sig2P <- invgamma::rinvgamma(S, shape = an/2, rate = dn/2)
PostRestCom <- mean(sapply(sig2P, function(x){dnorm(Brest, mean = bn[5], sd = (x*Bn[5,5])^0.5, log = FALSE)})) 
PostRest <- LaplacesDemon::dmvt(x = Brest, mu = bn[5], S = Hn[5,5], df = an, log = TRUE)
PostRest; log(PostRestCom)
# PostRestNew <- LaplacesDemon::dmvt(x = Brest, mu = bn[5], S = dn/an*Bn[5,5], df = an, log=TRUE)
# Prior evaluation
# Analytic
PriorRestNew <- LaplacesDemon::dmvt(x = Brest, mu = 0, S = cOpt*d0/a0, df = a0, log=TRUE)
PriorRestAn <- log(exp(lgamma((a0+1)/2))/((pi*cOpt*d0)^(1/2)*(exp(lgamma(a0/2)))))
# Computational
sig2 <- invgamma::rinvgamma(S, shape = a0/2, rate = d0/2)
PriorRestCom <- mean(sapply(sig2, function(x){dnorm(Brest, mean = 0, sd = (x*cOpt)^0.5, log = FALSE)})) 
PriorRestNew; PriorRestAn; log(PriorRestCom)
PriorRestAn/log(PriorRestCom)
PriorRestCom <- mean(sapply(sig2, function(x){dnorm(Brest, mean = 0, sd = (x*cOpt)^0.5, log = TRUE)})) 
BF12SD <- exp(PostRest - PriorRestNew)
2*log(1/BF12SD)
# Chib's method
sig2Post <- MCMCpack::rinvgamma(S,an/2,dn/2)
BetasGibbs <- sapply(1:S, function(s){MASS::mvrnorm(n = 1, mu = bn, Sigma = sig2Post[s]*Bn)})
# Mode function for continuous data
mode_continuous <- function(x){
  density_est <- density(x)               # Estimate density
  mode_value <- density_est$x[which.max(density_est$y)]  # Get x value at max density
  return(mode_value)
}
# Unrestricted model
BetasMode <- apply(BetasGibbs, 1, mode_continuous)
Sigma2Mode <- mode_continuous(sig2Post)
VarModel <- Sigma2Mode*diag(N)
MeanModel <- X%*%BetasMode
LogLik <- mvtnorm::dmvnorm(c(Y), mean = MeanModel, sigma = VarModel, log = TRUE, checkSymmetry = TRUE)
LogPrior <- mvtnorm::dmvnorm(BetasMode, mean = rep(0, K), sigma = Sigma2Mode*cOpt*diag(K), log = TRUE, checkSymmetry = TRUE)+log(MCMCpack::dinvgamma(Sigma2Mode, a0/2, d0/2))
LogPost1 <- mvtnorm::dmvnorm(BetasMode, mean = bn, sigma = Sigma2Mode*Bn, log = TRUE, checkSymmetry = TRUE)
LogPost2 <- log(MCMCpack::dinvgamma(Sigma2Mode, an/2, dn/2))
LogMarLikChib <- LogLik + LogPrior -(LogPost1 + LogPost2)
# Restricted model
anRest <- N + a0
XRest <- X[,-5]
KRest <- dim(XRest)[2]
B0Rest <- cOpt*diag(KRest) 
BnRest <- solve(solve(B0Rest)+t(XRest)%*%XRest)
bhatRest <- solve(t(XRest)%*%XRest)%*%t(XRest)%*%Y
b0Rest <- rep(0, KRest)
bnRest <- BnRest%*%(solve(B0Rest)%*%b0Rest+t(XRest)%*%XRest%*%bhatRest)
dnRest <- as.numeric(d0 + t(Y-XRest%*%bhatRest)%*%(Y-XRest%*%bhatRest)+t(bhatRest - b0Rest)%*%solve(solve(t(XRest)%*%XRest)+B0Rest)%*%(bhatRest - b0Rest))
sig2PostRest <- MCMCpack::rinvgamma(S,anRest/2,dnRest/2)
BetasGibbsRest <- sapply(1:S, function(s){MASS::mvrnorm(n = 1, mu = bnRest, Sigma = sig2PostRest[s]*BnRest)})
BetasModeRest <- apply(BetasGibbsRest, 1, mode_continuous)
Sigma2ModeRest <- mode_continuous(sig2PostRest)
VarModelRest <- Sigma2ModeRest*diag(N)
MeanModelRest <- XRest%*%BetasModeRest
LogLikRest <- mvtnorm::dmvnorm(c(Y), mean = MeanModelRest, sigma = VarModelRest, log = TRUE, checkSymmetry = TRUE)
LogPriorRest <- mvtnorm::dmvnorm(BetasModeRest, mean = rep(0, KRest), sigma = Sigma2ModeRest*cOpt*diag(KRest), log = TRUE, checkSymmetry = TRUE)+log(MCMCpack::dinvgamma(Sigma2ModeRest, a0/2, d0/2))
LogPost1Rest <- mvtnorm::dmvnorm(BetasModeRest, mean = bnRest, sigma = Sigma2ModeRest*BnRest, log = TRUE, checkSymmetry = TRUE)
LogPost2Rest <- log(MCMCpack::dinvgamma(Sigma2ModeRest, anRest/2, dnRest/2))
LogMarLikChibRest <- LogLikRest + LogPriorRest -(LogPost1Rest + LogPost2Rest)
BFChibs <- exp(LogMarLikChibRest-LogMarLikChib)
BFChibs; 1/BFChibs
2*log(1/BFChibs)
# Gelfand-Dey method
GDmarglik <- function(ids, X, Betas, MeanThetas, VarThetas, sig2Post){
  K <- dim(X)[2]
  Thetas <- c(Betas[ids,], sig2Post[ids])
  Lognom <- (1/(1-alpha))*mvtnorm::dmvnorm(Thetas, mean = MeanThetas, sigma = VarThetas, log = TRUE, checkSymmetry = TRUE)
  Logden1 <- mvtnorm::dmvnorm(Betas[ids,], mean = rep(0, K), sigma = sig2Post[ids]*cOpt*diag(K), log = TRUE, checkSymmetry = TRUE) + log(MCMCpack::dinvgamma(sig2Post[ids], a0/2, d0/2))
  VarModel <- sig2Post[ids]*diag(N)
  MeanModel <- X%*%Betas[ids,]
  Logden2 <- mvtnorm::dmvnorm(c(Y), mean = MeanModel, sigma = VarModel, log = TRUE, checkSymmetry = TRUE)
  LogGDid <- Lognom - Logden1 - Logden2
  return(LogGDid)
}
sig2Post <- MCMCpack::rinvgamma(S,an/2,dn/2)
Betas <- LaplacesDemon::rmvt(S, bn, Hn, an)
Thetas <- cbind(Betas, sig2Post)
MeanThetas <- colMeans(Thetas)
VarThetas <- var(Thetas)
iVarThetas <- solve(VarThetas)
ChiSQ <- sapply(1:S, function(s){(Thetas[s,]-MeanThetas)%*%iVarThetas%*%(Thetas[s,]-MeanThetas)})
alpha <- 0.01
criticalval <- qchisq(1-alpha, K + 1)
idGoodThetas <- which(ChiSQ <= criticalval)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
InvMargLik2 <- NULL
for(s in idGoodThetas){
  LogInvs <- GDmarglik(ids = s, X = X, Betas = Betas, MeanThetas = MeanThetas, VarThetas = VarThetas, sig2Post = sig2Post)
  InvMargLik2 <- c(InvMargLik2, LogInvs)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
summary(coda::mcmc(InvMargLik2))
mean(InvMargLik2)
# Restricted model
anRest <- N + a0
XRest <- X[,-5]
KRest <- dim(XRest)[2]
B0Rest <- cOpt*diag(KRest) 
BnRest <- solve(solve(B0Rest)+t(XRest)%*%XRest)
bhatRest <- solve(t(XRest)%*%XRest)%*%t(XRest)%*%Y
b0Rest <- rep(0, KRest)
bnRest <- BnRest%*%(solve(B0Rest)%*%b0Rest+t(XRest)%*%XRest%*%bhatRest)
dnRest <- as.numeric(d0 + t(Y-XRest%*%bhatRest)%*%(Y-XRest%*%bhatRest)+t(bhatRest - b0Rest)%*%solve(solve(t(XRest)%*%XRest)+B0Rest)%*%(bhatRest - b0Rest))
HnRest <- as.matrix(Matrix::forceSymmetric(dnRest*BnRest/anRest))
sig2PostRest <- MCMCpack::rinvgamma(S,anRest/2,dnRest/2)
BetasRest <- LaplacesDemon::rmvt(S, bnRest, HnRest, anRest)
ThetasRest <- cbind(BetasRest, sig2PostRest)
MeanThetasRest <- colMeans(ThetasRest)
VarThetasRest <- var(ThetasRest)
iVarThetasRest <- solve(VarThetasRest)
ChiSQRest <- sapply(1:S, function(s){(ThetasRest[s,]-MeanThetasRest)%*%iVarThetasRest%*%(ThetasRest[s,]-MeanThetasRest)})
idGoodThetasRest <- which(ChiSQRest <= criticalval)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
InvMargLik1 <- NULL
for(s in idGoodThetasRest){
  LogInvs <- GDmarglik(ids = s, X = XRest, Betas = BetasRest, MeanThetas = MeanThetasRest, VarThetas = VarThetasRest, sig2Post = sig2PostRest)
  InvMargLik1 <- c(InvMargLik1, LogInvs)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
summary(coda::mcmc(InvMargLik1))
mean(InvMargLik1)
BFFD <- exp(mean(InvMargLik2)-mean(InvMargLik1))
BFFD; mean(1/BFFD)
2*log(1/BFFD)
