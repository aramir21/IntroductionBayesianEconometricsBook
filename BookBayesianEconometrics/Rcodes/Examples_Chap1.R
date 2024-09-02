# Examples in chapter 1


########################## The Base rate fallacy ########################## 

PD <- 0.002 # Probability of disease
PPD <- 0.98 # True positive (Sensitivity)
PDP <- PD * PPD / (PD * PPD + (1 - PD) * (1 - PPD)) # Probability of disease given positive
paste("Probability of disease given a positive test is", sep = " ", round(PDP, 2))

##########################  Jesus Crist ########################## 

PR <- 2/(108.5 * 10^9) # Probability of resurrection
PWR <- 0.9999999 # Very reliable witness (true positive rate)
PRW <- PR * PWR / (PR * PWR + (1 - PR) * (1 - PWR)) # Probability of resurrection given witness
paste("Probability of resurrection given witness is", sep = " ", PRW)

##########################  The Monty Hall problem ########################## 
set.seed(0101) # Set simulation seed
S <- 100000 # Simulations
Game <- function(switch = 0){
  # switch = 0 is not change, and switch = 1 is to change
  opts <- 1:3 
  car <- sample(opts, 1) # car location
  guess1 <- sample(opts, 1) # Initial guess pick
  
  if(car != guess1) {
    host <- opts[-c(car, guess1)]
  } else {
    host <- sample(opts[-c(car, guess1)], 1)
  }
  
  win1 <- guess1 == car # Win given no change
  
  guess2 <- opts[-c(host, guess1)]
  
  win2 <- guess2 == car # Win given change
  
  if(switch == 0){
    win <- win1
  } else {
    win <- win2
  }
  
  return(win)
}

Prob <- mean(replicate(S, Game(switch = 0))) #Win probabilities not changing
paste("Winning probabilities no changing door is", Prob, sep = " ")

Prob <- mean(replicate(S, Game(switch = 1))) #Win probabilities changing
paste("Winning probabilities changing door is", Prob, sep = " ")

##########################  Health care insurance ########################## 

set.seed(010101)
y <- c(0, 3, 2, 1, 0) # Data
N <- length(y)

paste("The sample mean is", mean(y), sep = " ")

paste("The sample variance is", var(y), sep = " ")

ProbBo <- function(y, a0, b0){
  N <- length(y)
  an <- a0 + sum(y) # Posterior shape parameter
  bn <- b0 / ((b0 * N) + 1) # Posterior scale parameter
  p <- bn / (bn + 1) # Probability negative binomial density
  Pr <- 1 - pnbinom(0, size = an, prob = (1 - p)) # Probability of visiting the Doctor at least once next year
  # Observe that in R there is a slightly different parametrization.
  return(Pr)
} 

# Using a vague prior:

a0 <- 0.001 # Prior shape parameter
b0 <- 1 / 0.001 # Prior scale parameter
PriMeanV <- a0 * b0 # Prior mean
PriVarV <- a0 * b0^2 # Prior variance
paste("Prior mean and prior variance using vague information are", PriMeanV, "and", PriVarV, "respectively", sep = " ")

Pp <- ProbBo(y, a0 = 0.001, b0 = 1 / 0.001) # This setting is defining vague prior information.
paste("The probability of visiting the Doctor at least once next year using a vague prior is", Pp, sep = " ")

# Using Empirical Bayes
LogMgLik <- function(theta, y){
  N <- length(y) #sample size
  a0 <- theta[1] # prior shape hyperparameter
  b0 <- theta[2] # prior scale hyperparameter
  an <- sum(y) + a0 # posterior shape parameter
  if(a0 <= 0 || b0 <= 0){ #Avoiding negative values
    lnp <- -Inf
  }else{lnp <- lgamma(an) + sum(y)*log(b0/(N*b0+1)) - a0*log(N*b0+1) - lgamma(a0)} # log marginal likelihood
  return(-lnp)
}

theta0 <- c(0.01, 1/0.1) # Initial values
control <- list(maxit = 1000) # Number of iterations in optimization
EmpBay <- optim(theta0, LogMgLik, method = "BFGS", control = control, hessian = TRUE, y = y) # Optimization
EmpBay$convergence # Checking convergence

EmpBay$value # Maximum

a0EB <- EmpBay$par[1] # Prior shape using empirical Bayes
b0EB <- EmpBay$par[2] # Prior scale using empirical Bayes
paste("The prior shape and scale parameters are", a0EB, "and", b0EB, "respectively", sep = " ")

PriMeanEB <- a0EB * b0EB # Prior mean
PriVarEB <- a0EB * b0EB^2 # Prior variance
paste("Prior mean and variance using empirical Bayes are", PriMeanEB, "and", PriVarEB, "respectively", sep = " ")

PpEB <- ProbBo(y, a0 = a0EB, b0 = b0EB) # This setting is using emprical Bayes.
paste("The probability of visiting the Doctor at least once next year using empirical Bayes is", PpEB, sep = " ")

# Density figures
lambda <- seq(0.01, 10, 0.01) # Values of lambda
VaguePrior <- dgamma(lambda, shape = a0, scale = b0)
EBPrior <- dgamma(lambda, shape = a0EB, scale = b0EB)
PosteriorV <- dgamma(lambda, shape = a0 + sum(y), scale = b0 / ((b0 * N) + 1)) 
PosteriorEB <- dgamma(lambda, shape = a0EB + sum(y), scale = b0EB / ((b0EB * N) + 1))

# Likelihood function
Likelihood <- function(theta, y){
  LogL <- dpois(y, theta, log = TRUE)
  Lik <- prod(exp(LogL))
  return(Lik)
}
Liks <- sapply(lambda, function(par) {Likelihood(par, y = y)})
Sc <- max(PosteriorEB)/max(Liks) #Scale for displaying in figure
LiksScale <- Liks * Sc

data <- data.frame(cbind(lambda, VaguePrior, EBPrior, PosteriorV, PosteriorEB, LiksScale)) #Data frame

require(ggplot2) # Cool figures

require(latex2exp) # LaTeX equations in figures

require(ggpubr) # Multiple figures in one page

fig1 <- ggplot(data = data, aes(lambda, VaguePrior)) + 
  geom_line() +  
  xlab(TeX("$\\lambda$")) + ylab("Density") +
  ggtitle("Prior: Vague Gamma") 

fig2 <- ggplot(data = data, aes(lambda, EBPrior)) + 
  geom_line() +  
  xlab(TeX("$\\lambda$")) + ylab("Density") +
  ggtitle("Prior: Empirical Bayes Gamma")

fig3 <- ggplot(data = data, aes(lambda, PosteriorV)) + 
  geom_line() +  
  xlab(TeX("$\\lambda$")) + ylab("Density") +
  ggtitle("Posterior: Vague Gamma")

fig4 <- ggplot(data = data, aes(lambda, PosteriorEB)) + 
  geom_line() +  
  xlab(TeX("$\\lambda$")) + ylab("Density") +
  ggtitle("Posterior: Empirical Bayes Gamma")

FIG <- ggarrange(fig1, fig2, fig3, fig4,
                 ncol = 2, nrow = 2)

annotate_figure(FIG,
                top = text_grob("Vague versus Empirical Bayes: Poisson-Gamma model", color = "black", face = "bold", size = 14))

dataNew <- data.frame(cbind(rep(lambda, 3), c(EBPrior, PosteriorEB, LiksScale),
                            rep(1:3, each = 1000))) #Data frame

colnames(dataNew) <- c("Lambda", "Density", "Factor")
dataNew$Factor <- factor(dataNew$Factor, levels=c("1", "3", "2"), 
                         labels=c("Prior", "Likelihood", "Posterior"))

ggplot(data = dataNew, aes_string(x = "Lambda", y = "Density", group = "Factor")) + 
  geom_line(aes(color = Factor)) +
  xlab(TeX("$\\lambda$")) + ylab("Density") +
  ggtitle("Prior, likelihood and posterior: Empirical Bayes Poisson-Gamma model") +
  guides(color=guide_legend(title="Information")) +
  scale_color_manual(values = c("red", "yellow", "blue"))

# Predictive distributions
PredDen <- function(y, y0, a0, b0){
  N <- length(y)
  an <- a0 + sum(y) # Posterior shape parameter
  bn <- b0 / ((b0 * N) + 1) # Posterior scale parameter
  p <- bn / (bn + 1) # Probability negative binomial density
  Pr <- dnbinom(y0, size = an, prob = (1 - p)) # Predictive density
  # Observe that in R there is a slightly different parametrization.
  return(Pr)
}
y0 <- 0:10
PredVague <- PredDen(y = y, y0 = y0, a0 = a0, b0 = b0)
PredEB <- PredDen(y = y, y0 = y0, a0 = a0EB, b0 = b0EB)
dataPred <- as.data.frame(cbind(y0, PredVague, PredEB))
colnames(dataPred) <- c("y0", "PredictiveVague", "PredictiveEB")

ggplot(data = dataPred) + 
  geom_point(aes(y0, PredictiveVague, color = "red")) +  
  xlab(TeX("$y_0$")) + ylab("Density") +
  ggtitle("Predictive density: Vague and Empirical Bayes priors") +
  geom_point(aes(y0, PredictiveEB, color = "yellow")) +
  guides(color = guide_legend(title="Prior")) +
  scale_color_manual(labels = c("Vague", "Empirical Bayes"), values = c("red", "yellow")) +
  scale_x_continuous(breaks=seq(0,10,by=1))

# Posterior odds: Vague vs Empirical Bayes

PO12 <- exp(-LogMgLik(c(a0EB, b0EB), y = y))/exp(-LogMgLik(c(a0, b0), y = y))
paste("The posterior odds: Empirical Bayes vs Vague prior prior is", PO12, sep = " ")

PostProMEM <- PO12/(1 + PO12) # Posterior model probability Empirical Bayes
PostProbMV <- 1 - PostProMEM # Posterior model probability vague prior

paste("These are the posterior model probabilities", PostProMEM, PostProbMV, "for the Empirical Bayes and vague priors, respectively")

# Bayesian model average (BMA)
PostMeanEB <- (a0EB + sum(y)) * (b0EB / (b0EB * N + 1)) # Posterior mean Empirical Bayes 
PostMeanV <- (a0 + sum(y)) * (b0 / (b0 * N + 1)) # Posterior mean vague priors
BMAmean <- PostProMEM * PostMeanEB + PostProbMV * PostMeanV  # BMA posterior mean

PostVarEB <- (a0EB + sum(y)) * (b0EB / (b0EB * N + 1))^2 # Posterior variance Empirical Bayes
PostVarV <- (a0 + sum(y)) * (b0 / (b0 * N + 1))^2 # Posterior variance vague prior 

BMAVar <- PostProMEM * PostVarEB + PostProbMV * PostVarV + PostProMEM * (PostMeanEB - BMAmean)^2 + PostProbMV * (PostMeanV - BMAmean)^2# BMA posterior variance   

paste("The BMA posterior mean and variance are", BMAmean, "and", BMAVar, "respectively", sep = " ")

# BMA: Predictive
BMAPred <- PostProMEM * PredEB + PostProbMV * PredVague    
dataPredBMA <- as.data.frame(cbind(y0, BMAPred))
colnames(dataPredBMA) <- c("y0", "PredictiveBMA")

ggplot(data = dataPredBMA) + 
  geom_point(aes(y0, PredictiveBMA, color = "red")) +  
  xlab(TeX("$y_0$")) + ylab("Density") +
  ggtitle("Predictive density: BMA") +
  guides(color = guide_legend(title="BMA")) +
  scale_color_manual(labels = c("Probability"), values = c("red")) +
  scale_x_continuous(breaks=seq(0,10,by=1))

# Bayesian updating
BayUp <- function(y, lambda, a0, b0){
  N <- length(y)
  an <- a0 + sum(y) # Posterior shape parameter
  bn <- b0 / ((b0 * N) + 1) # Posterior scale parameter
  p <- dgamma(lambda, shape = an, scale = bn) # Posterior density
  return(list(Post = p, a0New = an, b0New = bn))
}

PostUp <- NULL
for(i in 1:N){
  if(i == 1){
    PostUpi <- BayUp(y[i], lambda, a0 = 0.001, b0 = 1/0.001)}
  else{
    PostUpi <- BayUp(y[i], lambda, a0 = PostUpi$a0New, b0 = PostUpi$b0New)
  }
  PostUp <- cbind(PostUp, PostUpi$Post)
}

DataUp <- data.frame(cbind(rep(lambda, 5), c(PostUp),
                           rep(1:5, each = 1000))) #Data frame

colnames(DataUp) <- c("Lambda", "Density", "Factor")

DataUp$Factor <- factor(DataUp$Factor, levels=c("1", "2", "3", "4", "5"), 
                        labels=c("Iter 1", "Iter 2", "Iter 3", "Iter 4", "Iter 5"))

ggplot(data = DataUp, aes_string(x = "Lambda", y = "Density", group = "Factor")) + 
  geom_line(aes(color = Factor)) +
  xlab(TeX("$\\lambda$")) + ylab("Density") +
  ggtitle("Bayesian updating: Poisson-Gamma model with vague prior") +
  guides(color=guide_legend(title="Update")) +
  scale_color_manual(values = c("red", "purple", "blue", "yellow", "black"))


S <- 100000 # Posterior draws
PostMeanLambdaUps <- sapply(1:N, function(i) {mean(sample(lambda, S, replace = TRUE, prob = PostUp[ , i]))}) #Posterior mean update i
paste("Posterior means using all information and sequential updating are:", round(PostMeanV, 2), "and", round(PostMeanLambdaUps[5], 2), sep = " ") 


PostVarLambdaUps <- sapply(1:N, function(i) {var(sample(lambda, S, replace = TRUE, prob = PostUp[ , i]))}) #Posterior variance update i
paste("Posterior variances using all information and sequential updating are:", round(PostVarV, 2), "and", round(PostVarLambdaUps[5], 2), sep = " ")


########################## Decision theory ########################## 

an <- sum(y) + a0EB # Posterior shape parameter
bn <- b0EB / (N*b0EB + 1) # Posterior scale parameter
S <- 1000000 # Number of posterior draws
Draws <- rgamma(1000000, shape = an, scale = bn) # Posterior draws
###### Point estimation ########
OptQua <- an*bn # Mean: Optimal choice quadratic loss function
OptQua

OptAbs <- qgamma(0.5, shape = an, scale = bn) # Median: Optimal choice absolute loss function
OptAbs

# Setting K0 = 2 and K1 = 1, that is, to underestimate lambda is twice as costly as to overestimate it.
K0 <- 2; K1 <- 1
OptGenAbs <- quantile(Draws, K0/(K0 + K1)) # Median: Optimal choice generalized absolute loss function
OptGenAbs

###### Hypothesis test ########
# H0: lambda in [0,1) vs H1: lambda in [1, Inf]
K0 <- 1; K1 <- 1
ProbH0 <- pgamma(1, shape = an, scale = bn) 
ProbH0 # Posterior  probability H0

ProbH1 <- 1 -ProbH0
ProbH1 # Posterior  probability H1

# we should reject H0 given ProbH1 > K1 / (K0 + K1) 

###### Credible intervals ########
LimInf <- qgamma(0.025, shape = an, scale = bn) # Lower bound
LimInf

LimSup <- qgamma(0.975, shape = an, scale = bn) # Upper bound
LimSup

HDI <- HDInterval::hdi(Draws, credMass = 0.95) # Highest posterior density credible interval
HDI

###### Predictive optimal choices ########
p <- bn / (bn + 1) # Probability negative binomial density
OptPred <- p/(1-p)*an # Optimal point prediction given a quadratic loss function in prediction
OptPred

