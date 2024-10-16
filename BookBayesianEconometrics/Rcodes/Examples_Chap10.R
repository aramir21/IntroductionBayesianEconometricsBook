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
df <- as.data.frame(cbind(y, X)) 
colnames(df) <- c("y", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
write.csv(df, file="51SimNormalBMANew.csv", row.names=FALSE)
#### BIC approximation
BMAglm <- BMA::bic.glm(y ~ X, data = df, glm.family = gaussian()) 
summary(BMAglm)

#### Markov chain Monte Carlo model composition
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
#### From scratch
combs <- expand.grid(c(0,1), c(0,1), c(0,1), c(0,1), c(0,1),c(0,1), c(0,1), c(0,1), c(0,1), c(0,1))
mll <- rep(NA, 2^K)
Xnew <- apply(X, 2, scale)
vpos <- N
mll <- rep(NA, 2^K)
# Initial exact exercise
for (i in 1:(2^K)) {
  indr <- combs[i, ] == 1
  kr <- sum(indr)
  if(kr == 0){
    gr <- ifelse(N > kr^2, 1/N, kr^(-2))
    PX <- diag(N)
    s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
    mll[i] <- (kr/2)*log(gr/(1+gr))-(N-1)/2*log(s2pos)
  }else{
    gr <- ifelse(N > kr^2, 1/N, kr^(-2))
    Xr <- matrix(Xnew[ , indr], ncol = kr)
    PX <- diag(N) - Xr%*%solve(t(Xr)%*%Xr)%*%t(Xr)
    s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
    mll[i] <- (kr/2)*log(gr/(1+gr))-(N-1)/2*log(s2pos)
  }
}
# Results
MaxPMP <- which.max(mll)
StMarLik <- exp(mll-max(mll))
PMP <- StMarLik/sum(StMarLik) 
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

# MC3 methods
# Initial parameters
final <- 500
burnin <- 0.2*final
iter <- burnin + final

# M best models
M <- 1000
vals <- matrix(runif(K*M), ncol = K, nrow = M)
model.mat <- ifelse(vals <= 0.5, 0, 1)
if(any(rowSums(model.mat == rep(0, K)) == 5)) {
  model.mat[which(rowSums(model.mat == rep(0, K)) == 5), ] <- rep(1, K)
}
mmll <- rep(NA, M)
for (i in 1:M) {
  ind <- (model.mat[i, ] == 1)
  kr <- sum(ind)
  if(kr == 0){
    gr <- ifelse(N > kr^2, 1/N, kr^(-2))
    PX <- diag(N)
    s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
    mmll[i] <- (kr/2)*log(gr/(1+gr))-(N-1)/2*log(s2pos)
  }else{
    gr <- ifelse(N > kr^2, 1/N, kr^(-2))
    Xr <- matrix(Xnew[ , ind], ncol = kr)
    PX <- diag(N) - Xr%*%solve(t(Xr)%*%Xr)%*%t(Xr)
    s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
    mmll[i] <- (kr/2)*log(gr/(1+gr))-(N-1)/2*log(s2pos)
  }
}
StMarLik <- exp(mmll-max(mmll))
nmmll <- StMarLik/sum(StMarLik) 
oind <- order(nmmll, decreasing = TRUE)
mmll <- mmll[oind]
model.mat <- model.mat[oind, ]

for (i in 1:iter) {
  indi <- (model.mat[M, ] == 1)
  kr <- sum(indi)
  repeat{
    vals <- runif(K)
    cmodel <- ifelse(vals <= 0.5, 0, 1)
    indc <- (cmodel == 1)
    if (abs(kr - sum(indc)) <= 1 & sum(indc) != 0) {
      break
    }
  }
  Xr <- matrix(Xnew[ , indi], ncol = kr)
  PX <- diag(N) - Xr%*%solve(t(Xr)%*%Xr)%*%t(Xr)
  s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
  imll <- (kr/2)*log(gr/(1+gr))-(N-1)/2*log(s2pos)
  
  Xr <- matrix(Xnew[ , indc], ncol = sum(indc))
  PX <- diag(N) - Xr%*%solve(t(Xr)%*%Xr)%*%t(Xr)
  s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
  cmll <- (kr/2)*log(gr/(1+gr))-(N-1)/2*log(s2pos)
  
  alpha <- min(exp(cmll-imll), 1)
  U <- runif(1)
  if (U <= alpha) {
    model.mat[M, ] <- cmodel
    mmll[M] <- cmll
  }
  StMarLik <- exp(mmll-max(mmll))
  nmmll <- StMarLik/sum(StMarLik) 
  oind <- order(nmmll, decreasing = TRUE)
  mmll <- mmll[oind]
  model.mat <- model.mat[oind, ]
}

(pip <- colMeans(model.mat))

# Initial model
imodel <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
final.models <- matrix(NA, nrow = iter, ncol = K)

for (i in 1:iter) {
  indi <- (imodel == 1)
  kr <- sum(indi)
  repeat{
    vals <- runif(K)
    cmodel <- ifelse(vals <= 0.5, 0, 1)
    indc <- (cmodel == 1)
    if (abs(kr - sum(indc)) <= 1 & sum(indc) != 0) {
      break
    }
  }
  Xr <- matrix(Xnew[ , indi], ncol = kr)
  PX <- diag(N) - Xr%*%solve(t(Xr)%*%Xr)%*%t(Xr)
  s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
  imll <- log((gr/(1+gr))^(kr/2))*(s2pos^(-(N-1)/2))
  
  Xr <- matrix(Xnew[ , indc], ncol = sum(indc))
  PX <- diag(N) - Xr%*%solve(t(Xr)%*%Xr)%*%t(Xr)
  s2pos <- c(t(y)%*%PX%*%y/(1 + gr) + gr*(t(y - mean(y))%*%(y - mean(y)))/(1 + gr))
  cmll <- log((gr/(1+gr))^(sum(indc)/2))*(s2pos^(-(N-1)/2))
  
  alpha <- min(exp(cmll-imll), 1)
  U <- runif(1)
  if (U <= alpha) {
    imodel <- cmodel
  }
  final.models[i, ] <- imodel
}

final.models <- final.models[-(1:burnin), ]
(pip <- colMeans(final.models))



