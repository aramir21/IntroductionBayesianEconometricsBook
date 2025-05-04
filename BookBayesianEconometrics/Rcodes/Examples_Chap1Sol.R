# Examples in chapter 1

######################### Cab companies ######################################
B <- 0.15 # Blue cabs
G <- 0.85 # Green cabs
TP <- 0.8 # True positive witness
FP <- 0.25 # False positive
BW <- TP * B / (TP * B + FP * (1 - B))
paste("Probability of a blue cab given witness is", sep = " ", round(BW, 2))


##########################  The Monty Hall problem reload ########################## 
set.seed(0101) # Set simulation seed
S <- 100000 # Simulations
Game <- function(opt = 3){
  # opt: number of options. opt > 2
  opts <- 1:opt 
  car <- sample(opts, 1) # car location
  guess1 <- sample(opts, 1) # Initial guess
  
  if(opt == 3 && car != guess1) {
    host <- opts[-c(car, guess1)]
  } else {
    host <- sample(opts[-c(car, guess1)], 1)
  }
  
  win1 <- guess1 == car # Win given no change
  
  if(opt == 3) {
    guess2 <- opts[-c(host, guess1)]
  } else {
    guess2 <- sample(opts[-c(host, guess1)], 1)
  } 
  
  win2 <- guess2 == car # Win given change
  
  return(c(win1, win2))
}
#Win probabilities
Prob <- rowMeans(replicate(S, Game(opt = 4))) 
#Winning probabilities no changing door
Prob[1]
#Winning probabilities changing door
Prob[2]

##########################  Health insurance reload ##########################
set.seed(010101)
y <- c(0, 3, 2, 1, 0) # Data
N <- length(y)

# Predictive distribution
ProbBo <- function(y, a0, b0){
  N <- length(y) 
  #sample size
  aN <- a0 + sum(y) 
  # Posterior shape parameter
  bN <- b0 + N 
  # Posterior scale parameter
  p <- 1 / (bN + 1) 
  # Probability negative binomial density
  Pr <- 1 - pnbinom(0, size = aN, prob = (1 - p)) 
  # Probability of visiting the Doctor
  # Observe that in R there is a slightly 
  # different parametrization.
  return(Pr)
} 

# Using a vague prior:
a0 <- 0.001 # Prior shape parameter
b0 <- 0.001 # Prior scale parameter
PriMeanV <- a0 / b0 # Prior mean
PriVarV <- a0 / b0^2 # Prior variance
Pp <- ProbBo(y, a0 = 0.001, b0 = 0.001) 
# This setting is vague prior information.
Pp

# Using Emprirical Bayes
LogMgLik <- function(theta, y){
  N <- length(y) 
  #sample size
  a0 <- theta[1] 
  # prior shape hyperparameter
  b0 <- theta[2] 
  # prior scale hyperparameter
  aN <- sum(y) + a0 
  # posterior shape parameter
  if(a0 <= 0 || b0 <= 0){ 
    #Avoiding negative values
    lnp <- -Inf
  }else{lnp <- lgamma(aN) - sum(y)*log(b0+N) + 
    a0*log(b0/(b0+N)) - lgamma(a0)} 
  # log marginal likelihood
  return(-lnp)
}

theta0 <- c(0.01, 0.01) 
# Initial values
control <- list(maxit = 1000) 
# Number of iterations in optimization
EmpBay <- optim(theta0, LogMgLik, method = "BFGS", 
                control = control, hessian = TRUE, y = y) 
# Optimization
EmpBay$convergence 
# Checking convergence
EmpBay$value # Maximum
a0EB <- EmpBay$par[1] 
# Prior shape using empirical Bayes
a0EB
b0EB <- EmpBay$par[2] 
# Prior scale using empirical Bayes
b0EB 

PriMeanEB <- a0EB / b0EB 
# Prior mean
PriVarEB <- a0EB / b0EB^2 
# Prior variance
PpEB <- ProbBo(y, a0 = a0EB, b0 = b0EB) 
# This setting is using empirical Bayes.
PpEB

##########################  Car claim ##########################
set.seed(010101)
y <- c(0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0) 
# Data
N <- length(y)

#require(TailRank)

# Predictive distribution
ProbBo <- function(y, a0, b0){
  N <- length(y) 
  #sample size
  aN <- a0 + sum(y) 
  # Posterior shape parameter
  bN <- b0 + N - sum(y) 
  # Posterior scale parameter
  pr <- aN / (aN + bN) 
  # Probability of a claim the next year
  return(pr)
} 

# Using a vague prior:
a0 <- 1 # Prior shape parameter
b0 <- 1 # Prior scale parameter
PriMeanV <- a0 / (a0 + b0) 
# Prior mean
PriVarV <- (a0*b0) / (((a0+b0)^2)*(a0+b0+1)) 
# Prior variance

Pp <- ProbBo(y, a0 = 1, b0 = 1)
# This setting is defining vague prior information.
# The probability of a claim 
Pp
# Using Emprirical Bayes
LogMgLik <- function(theta, y){
  N <- length(y) 
  #sample size
  a0 <- theta[1] 
  # prior shape hyperparameter
  b0 <- theta[2] 
  # prior scale hyperparameter
  aN <- sum(y) + a0 
  # posterior shape parameter
  if(a0 <= 0 || b0 <= 0){ 
    #Avoiding negative values
    lnp <- -Inf
  }else{lnp <- lgamma(a0+b0) + lgamma(aN) + 
    lgamma(b0+N-sum(y)) -lgamma(a0) -lgamma(b0) - 
    lgamma(a0+b0+N)} 
  # log marginal likelihood
  return(-lnp)
}

theta0 <- c(0.1, 0.1) 
# Initial values
control <- list(maxit = 1000) 
# Number of iterations in optimization
EmpBay <- optim(theta0, LogMgLik, method = "BFGS", 
                control = control, hessian = TRUE, y = y) 
# Optimization
EmpBay$convergence 
# Checking convergence
EmpBay$value # Maximum
a0EB <- EmpBay$par[1] 
# Prior shape using empirical Bayes
b0EB <- EmpBay$par[2] 
# Prior scale using empirical Bayes

PriMeanEB <- a0EB /(a0EB + b0EB)  
# Prior mean
PriVarEB <- (a0EB*b0EB)/(((a0EB+b0EB)^2)*(a0EB+b0EB+1)) 
# Prior variance

PpEB <- ProbBo(y, a0 = a0EB, b0 = b0EB)
# This setting is using empirical Bayes.
PpEB
# Density figures
lambda <- seq(0.001, 1, 0.001) 
# Values of lambda
VaguePrior <- dbeta(lambda, shape1 = a0, shape2 = b0)
EBPrior <- dbeta(lambda, shape1 = a0EB, shape2 = b0EB)
PosteriorV <- dbeta(lambda, shape1 = a0 + sum(y), 
                    shape2 = b0 + N - sum(y))
PosteriorEB <- dbeta(lambda, shape1 = a0EB + sum(y),
                     shape2 = b0EB + N - sum(y))

# Likelihood function
Likelihood <- function(theta, y){
  LogL <- dbinom(y, 1, theta, log = TRUE)
  #  LogL <- dbern(y, theta)
  Lik <- prod(exp(LogL))
  return(Lik)
}


Liks <- sapply(lambda, function(par) {
  Likelihood(par, y = y)})
Sc <- max(PosteriorEB)/max(Liks) 
#Scale for displaying in figure
LiksScale <- Liks * Sc

data <- data.frame(cbind(lambda, VaguePrior, EBPrior, 
                         PosteriorV, PosteriorEB, LiksScale)) 
#Data frame

require(ggplot2) 
# Cool figures
require(latex2exp) 
# LaTeX equations in figures
require(ggpubr) 
# Multiple figures in one page

fig1 <- ggplot(data = data, aes(lambda, VaguePrior)) + 
  geom_line() +  
  xlab(TeX("$p$")) + ylab("Density") + 
  ggtitle("Prior: Vague Beta") 

fig2 <- ggplot(data = data, aes(lambda, EBPrior)) + 
  geom_line() +  
  xlab(TeX("$p$")) + ylab("Density") +
  ggtitle("Prior: Empirical Bayes Beta")

fig3 <- ggplot(data = data, aes(lambda, PosteriorV)) + 
  geom_line() +  
  xlab(TeX("$p$")) + ylab("Density") +
  ggtitle("Posterior: Vague Beta")

fig4 <- ggplot(data = data, aes(lambda, PosteriorEB)) + 
  geom_line() +  
  xlab(TeX("$p$")) + ylab("Density") +
  ggtitle("Posterior: Empirical Bayes Beta")

FIG <- ggarrange(fig1, fig2, fig3, fig4,
                 ncol = 2, nrow = 2)

annotate_figure(FIG,
                top = text_grob("Vague versus Empirical Bayes: 
Beta-Bernoulli model", color = "black", 
                                face = "bold", size = 14))

# Prior, likelihood and posterior: 
#Empirical Bayes Binonial-Beta model
dataNew <- data.frame(cbind(rep(lambda, 3), 
                            c(EBPrior, PosteriorEB, LiksScale),
                            rep(1:3, each = 1000))) 
#Data frame

colnames(dataNew) <- c("Lambda", "Density", "Factor")
dataNew$Factor <- factor(dataNew$Factor, levels=c("1", "3", 
                                                  "2"), labels=c("Prior", "Likelihood", "Posterior"))

ggplot(data = dataNew, aes_string(x = "Lambda", 
                                  y = "Density", group = "Factor")) + 
  geom_line(aes(color = Factor)) +
  xlab(TeX("$\\lambda$")) + ylab("Density") +
  ggtitle("Prior, likelihood and posterior: Empirical Bayes
 Poisson-Gamma model") +
  guides(color=guide_legend(title="Information")) +
  scale_color_manual(values = c("red", "yellow", "blue"))

# Predictive distributions
require(TailRank)

PredDen <- function(y, y0, a0, b0){
  N <- length(y)
  aN <- a0 + sum(y) # Posterior shape parameter
  bN <- b0 + N - sum(y) # Posterior scale parameter
  Pr <- aN/(aN+bN)
  Probs <- dbinom(y0, 1, prob = Pr)
  return(Probs)
}
y0 <- 0:1
PredVague <- PredDen(y = y, y0 = y0, a0 = a0, b0 = b0)
PredEB <- PredDen(y = y, y0 = y0, a0 = a0EB, b0 = b0EB)
dataPred <- as.data.frame(cbind(y0, PredVague, PredEB))
colnames(dataPred) <- c("y0", "PredictiveVague", 
                        "PredictiveEB")

ggplot(data = dataPred) + 
  geom_point(aes(y0, PredictiveVague, color = "red")) +  
  xlab(TeX("$y_0$")) + ylab("Density") +
  ggtitle("Predictive density: Vague and Empirical Bayes priors") + 
  geom_point(aes(y0, PredictiveEB, color = "yellow")) +
  guides(color = guide_legend(title="Prior")) +
  scale_color_manual(labels = c("Vague", "Empirical Bayes"), 
                     values = c("red", "yellow")) +
  scale_x_continuous(breaks=seq(0,1,by=1))

# Posterior odds: Vague vs Empirical Bayes
PO12 <- exp(-LogMgLik(c(a0EB, b0EB), y = y))/
  exp(-LogMgLik(c(a0, b0), y = y))

PostProMEM <- PO12/(1 + PO12) 
# Posterior model probability Empirical Bayes
PostProMEM
PostProbMV <- 1 - PostProMEM 
# Posterior model probability vague prior
PostProbMV

# Bayesian model average (BMA)
PostMeanEB <- (a0EB + sum(y)) / (a0EB + b0EB + N) 
# Posterior mean Empirical Bayes 
PostMeanV <- (a0 + sum(y)) / (a0 + b0 + N) 
# Posterior mean vague priors
BMAmean <- PostProMEM * PostMeanEB + PostProbMV * PostMeanV  
# BMA posterior mean

PostVarEB <- (a0EB + sum(y))*(b0EB + N - sum(y)) / 
  ((a0EB + b0EB + N)^2)*(a0EB + b0EB + N -1) 
# Posterior variance Empirical Bayes
PostVarV <- (a0 + sum(y))*(b0 + N - sum(y)) / 
  ((a0 + b0 + N)^2)*(a0 + b0 + N -1) 
# Posterior variance vague prior 

BMAVar <- PostProMEM * PostVarEB + PostProbMV * PostVarV 
+ PostProMEM * (PostMeanEB - BMAmean)^2 + PostProbMV * 
  (PostMeanV - BMAmean)^2
# BMA posterior variance   

# BMA: Predictive
BMAPred <- PostProMEM * PredEB + PostProbMV * PredVague    
dataPredBMA <- as.data.frame(cbind(y0, BMAPred))
colnames(dataPredBMA) <- c("y0", "PredictiveBMA")

ggplot(data = dataPredBMA) + 
  geom_point(aes(y0, PredictiveBMA, color = "red")) +  
  xlab(TeX("$y_0$")) + ylab("Density") +
  ggtitle("Predictive density: BMA") +
  guides(color = guide_legend(title="BMA")) +
  scale_color_manual(labels = c("Probability"), 
                     values = c("red")) + 
  scale_x_continuous(breaks=seq(0,1,by=1))

# Bayesian updating
BayUp <- function(y, lambda, a0, b0){
  N <- length(y)
  aN <- a0 + sum(y) 
  # Posterior shape parameter
  bN <- b0 + N - sum(y)    
  # Posterior scale parameter
  p <- dbeta(lambda, shape1 = aN, shape2 = bN) 
  # Posterior density
  return(list(Post = p, a0New = aN, b0New = bN))
}

PostUp <- NULL
for(i in 1:N){
  if(i == 1){
    PostUpi <- BayUp(y[i], lambda, a0 = 1, b0 = 1)}
  else{
    PostUpi <- BayUp(y[i], lambda, 
                     a0 = PostUpi$a0New, b0 = PostUpi$b0New)
  }
  PostUp <- cbind(PostUp, PostUpi$Post)
}

DataUp <- data.frame(cbind(rep(lambda, 15), c(PostUp),
                           rep(1:15, each = 1000))) 
#Data frame
colnames(DataUp) <- c("Lambda", "Density", "Factor")
DataUp$Factor <- factor(DataUp$Factor, levels=c("1","2",
                                                "3","4","5","6","7","8","9","10","11","12","13","14","15"), 
                        labels=c("Iter_1","Iter_2","Iter_3","Iter_4","Iter_5",
                                 "Iter_6","Iter_7","Iter_8","Iter_9","Iter_10","Iter_11",
                                 "Iter_12","Iter_13","Iter_14","Iter_15"))

ggplot(data = DataUp, aes_string(x = "Lambda", 
                                 y = "Density", group = "Factor")) + 
  geom_line(aes(color = Factor)) +
  xlab(TeX("$p$")) + ylab("Density") +
  ggtitle("Bayesian updating: 
 Beta-Binomial model with vague prior") +
  guides(color=guide_legend(title="Update")) 		







