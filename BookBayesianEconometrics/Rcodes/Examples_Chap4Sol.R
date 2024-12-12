rm(list = ls())
set.seed(010101)
##########################  Multinomial-Dirichlet exercise: Polling 2016 USA presidential race ########################## 
y <- c(411, 373, 149) 
# Clinton, Trump, Other
# Public Policy Polling September 27-28, 2016 five-way race
# https://en.wikipedia.org/wiki/Nationwide_opinion_polling_for_the_2016_United_States_presidential_election#Aggregate_polls
alpha0 <- rep(1, 3)
# Hyperparameters: non-informative distribution
alphan <- alpha0 + y 
S <- 100000 
# Sample draws
thetas <- MCMCpack::rdirichlet(S, alphan)
colnames(thetas) <- c("Clinton", "Trump", "Other")
head(thetas)
# Compositional::bivt.contour(thetas, cont.line = FALSE, appear = FALSE)
dif <- thetas[,1] - thetas[,2]
data <- data.frame(dif)
names(data) <- c("Difference")
library(ggplot2)
p <- ggplot(data) +  
  geom_histogram(aes(x = Difference), binwidth = 0.01) +
  geom_vline(xintercept=0.0, lwd=1, colour="red") + 
  ggtitle("Percentage difference Clinton vs Trump 2016 presidential race") +
  xlab("Percentage Difference") + ylab("")
p
difmcmc <- coda::mcmc(dif)
summary(difmcmc)
CW <- mean(difmcmc>0)

# Predictive
y0 <- round(alphan/sum(alphan)*100,0) # c(44, 40, 16)
Pred <- apply(thetas, 1, function(p) {rmultinom(1, size = sum(y0), prob = p)})
sum(sapply(1:S, function(s) {sum(Pred[,s] == y0) == 3}))/S

PredY0 <- function(y0){
  n <- sum(y0)
  Res1 <- sum(sapply(1:length(y), function(l){lgamma(alphan[l] + y0[l]) - lgamma(alphan[l]) - lfactorial(y0[l])}))
  Res <- lfactorial(n) + lgamma(sum(alphan)) - lgamma(sum(alphan)+n) + Res1
  return(exp(Res))
}
PredY0(y0)

########################## Math test example ########################## 
rm(list =)
set.seed(010101)
N <- 50
# Sample size
muhat <- 102
# Sample mean
sig2hat <- 10
# Sample variance

# Hyperparameters
mu0 <- 100
beta0 <- 1
delta0 <- 0.001
alpha0 <- 0.001

S <- 100000
# Posterior draws
alphan <- alpha0 + N
deltan <- sig2hat*(N - 1) + delta0 + beta0*N/(beta0 + N)*(muhat - mu0)^2
sig2Post <- invgamma::rinvgamma(S, shape = alphan, rate = deltan)
summary(sig2Post)
betan <- beta0 + N
mun <- (beta0*mu0 + N*muhat)/betan
# muPost <- sapply(sig2Post, function(s2){rnorm(1, mun, sd = (s2/betan)^0.5)})
# 
# muPostq <- quantile(muPost, c(0.025, 0.5, 0.975))
# muPostq
# cutoff <- 103
# PmuPostcutoff <- mean(muPost > cutoff)
# PmuPostcutoff

# Using Student's t
muPost_t <- ((deltan/(alphan*betan))^0.5)*rt(S, alphan) + mun
c1 <- rgb(173,216,230,max = 255, alpha = 50, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 50, names = "lt.pink")
# hist(muPost, main = "Histogram: Posterior mean", xlab = "Posterior mean", col = c2)
# hist(muPost_t, main = "Histogram: Posterior mean", xlab = "Posterior mean", add = T, col = c1)
hist(muPost_t, main = "Histogram: Posterior mean", xlab = "Posterior mean", col = c1)
muPost_tq <- quantile(muPost_t, c(0.025, 0.5, 0.975))
muPost_tq
PmuPost_tcutoff <- mean(muPost_t > cutoff)
PmuPost_tcutoff

########################## Empirical Bayes: Linear regression (electricity demand) ########################## 
rm(list = ls())
set.seed(010101)
# Electricity demand
DataUt <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/Utilities.csv", sep = ",", header = TRUE, quote = "")
DataUtEst <- DataUt %>% 
  filter(Electricity != 0)
attach(DataUtEst)
# Dependent variable: Monthly consumption (kWh) in log
Y <- log(Electricity)
N <- length(Y)
# Regressors quantity including intercept
X <- cbind(LnPriceElect, IndSocio1, IndSocio2, Altitude, Nrooms, HouseholdMem, Children, Lnincome, 1)
# Regressor without price
Xnew <- cbind(IndSocio1, IndSocio2, Altitude, Nrooms, HouseholdMem, Children, Lnincome, 1)
# Log marginal function (multiply by -1 due to minimization)
LogMarLikLM <- function(X, c0){
  k <- dim(X)[2]
  N <- dim(X)[1]	
  # Hyperparameters
  B0 <- c0*diag(k)
  b0 <- rep(0, k)
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
# Hyperparameters
d0 <- 0.001
a0 <- 0.001
# Empirical Bayes: Obtain c0 maximizing the log marginal likelihood
c0 <- 1000 
EB <- optim(c0, fn = LogMarLikLM, method = "Brent", lower = 0.0001, upper = 10^6, X = X)
EB$par
EB$value
EBnew <- optim(c0, fn = LogMarLikLM, method = "Brent", lower = 0.0001, upper = 10^6, X = Xnew)
EBnew$par
EBnew$value
# Change of order to take into account the -1 in the LogMarLikLM function
BFEM <- exp(EBnew$value - EB$value) 
BFEM

########################## Utilities demand: Multivariate regression ########################## 
rm(list = ls())
set.seed(010101)
library(dplyr)
# Electricity demand
DataUt <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/Utilities.csv", sep = ",", header = TRUE, quote = "")
DataUtEst <- DataUt %>%  
  filter(Electricity != 0 & Water !=0 & Gas != 0)
attach(DataUtEst)
Y <- cbind(log(Electricity), log(Water), log(Gas))
X <- cbind(LnPriceElect, LnPriceWater, LnPriceGas, IndSocio1, IndSocio2, Altitude, Nrooms, HouseholdMem, Children, Lnincome, 1)
M <- dim(Y)[2]
K <- dim(X)[2]
N <- dim(Y)[1]
# Hyperparameters
B0 <- matrix(0, K, M)
c0 <- 1000
V0 <- c0*diag(K)
Psi0 <- c0*diag(M)
a0 <- M
# Posterior parameters
Bhat <- solve(t(X)%*%X)%*%t(X)%*%Y 
S <- t(Y - X%*%Bhat)%*%(Y - X%*%Bhat)
Vn <- solve(solve(V0) + t(X)%*%X) 
Bn <- Vn%*%(solve(V0)%*%B0 + t(X)%*%X%*%Bhat)
Psin <- Psi0 + S + t(B0)%*%solve(V0)%*%B0 + t(Bhat)%*%t(X)%*%X%*%Bhat - t(Bn)%*%solve(Vn)%*%Bn
an <- a0 + N
#Posterior draws
s <- 10000 #Number of posterior draws
SIGs <- replicate(s, LaplacesDemon::rinvwishart(an, Psin))
# BsCond <- sapply(1:s, function(s) {MixMatrix::rmatrixnorm(n = 1, mean=Bn, U = Vn,V = SIGs[,,s])})
# summary(coda::mcmc(t(BsCond)))
Bs <- sapply(1:s, function(s) {MixMatrix::rmatrixt(n = 1, mean=Bn, U = Vn,V = Psin, df = an + 1 - M)})
summary(coda::mcmc(t(Bs)))
SIGMs <- t(sapply(1:s, function(l) {gdata::lowerTriangle(SIGs[,,l], diag=TRUE, byrow=FALSE)}))
summary(coda::mcmc(SIGMs))
hdiBs <- HDInterval::hdi(t(Bs), credMass = 0.95) # Highest posterior density credible interval
hdiBs
hdiSIG <- HDInterval::hdi(SIGMs, credMass = 0.95) # Highest posterior density credible interval
hdiSIG
# Log marginal function (multiply by -1 due to minimization)
LogMarLikLM <- function(X, c0){
  c10 <- c0[1]; c20 <- c0[2]
  k <- dim(X)[2]
  N <- dim(X)[1]
  # Hyperparameters
  V0 <- c10*diag(K)
  Psi0 <- c20*diag(M)
  # Posterior parameters
  Bhat <- solve(t(X)%*%X)%*%t(X)%*%Y 
  S <- t(Y - X%*%Bhat)%*%(Y - X%*%Bhat)
  Vn <- solve(solve(V0) + t(X)%*%X) 
  Bn <- Vn%*%(solve(V0)%*%B0 + t(X)%*%X%*%Bhat)
  Psin <- Psi0 + S + t(B0)%*%solve(V0)%*%B0 + t(Bhat)%*%t(X)%*%X%*%Bhat - t(Bn)%*%solve(Vn)%*%Bn
  # Log marginal likelihood
  logpy <- (N*M/2)*log(1/pi)+(a0/2)*log(det(Psi0)) - (an/2)*log(det(Psin)) + (M/2)*(log(det(Vn)) - log(det(V0))) + lgamma(an/2)-lgamma(a0/2)
  return(-logpy)
}
c0 <- rep(1000, 2)
LogML <- LogMarLikLM(X=X, c0 = c0)
# Using income in dollars as regressor
Xnew <- cbind(LnPriceElect, LnPriceWater, LnPriceGas, IndSocio1, IndSocio2, Altitude, Nrooms, HouseholdMem, Children, exp(Lnincome), 1)
LogMLnew <- LogMarLikLM(X=Xnew, c0 = c0)
# Bayes factor
BF12 <- exp(LogMLnew - LogML)
BF12
# Using income in thousand dollars as regressor
XnewT <- cbind(LnPriceElect, LnPriceWater, LnPriceGas, IndSocio1, IndSocio2, Altitude, Nrooms, HouseholdMem, Children, exp(Lnincome)/1000, 1)
LogMLnewT <- LogMarLikLM(X=XnewT, c0 = c0)
# Bayes factor
BF13 <- exp(LogMLnewT - LogML)
BF13
# Empirical Bayes: Obtain c0 maximizing the log marginal likelihood
EB <- optim(c0, fn = LogMarLikLM, method = "BFGS", X = X)
EB$par
EB$value
EBnew <- optim(c0, fn = LogMarLikLM, method = "BFGS", X = Xnew)
EBnew$par
EBnew$value
# Bayes factor
BF12EB <- exp(EBnew$value - EB$value)
BF12EB
# Empirical Bayes: Obtain c0 maximizing the log marginal likelihood
EBnewT <- optim(c0, fn = LogMarLikLM, method = "BFGS", X = XnewT)
EBnewT$par
EBnewT$value
# Bayes factor
BF13EB <- exp(EBnewT$value - EB$value)
BF13EB
# Predictive distribution
Xpred <- c(log(0.15), log(0.70), log(0.75), 1, 0, 0, 2, 3, 1, log(500), 1)
Mean <- Xpred%*%Bn
Hn <- 1+t(Xpred)%*%Vn%*%Xpred
UtilDemand <- exp(replicate(s, MixMatrix::rmatrixt(n = 1, mean=Mean, U = Hn, V = Psin, df = an + 1 - M)))
ElePred <- UtilDemand[1,1,]
WatPred <- UtilDemand[1,2,]
GasPred <- UtilDemand[1,3,]

data <- data.frame(cbind(ElePred, WatPred, GasPred)) #Data frame
annotations1 <- data.frame(
  x = round(quantile(data$ElePred, c(0.025, 0.5, 0.975)),1),
  y = c(600, 1000, 600),
  label = c("2.5%:", "50%:", "97.5%:")
)

annotations2 <- data.frame(
  x = round(quantile(data$WatPred, c(0.025, 0.5, 0.975)),1),
  y = c(600, 1000, 600),
  label = c("2.5%:", "50%:", "97.5%:")
)

annotations3 <- data.frame(
  x = round(quantile(data$GasPred, c(0.025, 0.5, 0.975)),1),
  y = c(600, 1000, 600),
  label = c("2.5%:", "50%:", "97.5%:")
)

require(ggplot2) # Cool figures
require(ggpubr) # Multiple figures in one page
require(latex2exp) # LaTeX equations in figures
fig1 <- ggplot(data = data, aes(ElePred)) + geom_histogram(bins = 40, color = "#000000", fill = "#0099F8") + 	xlab("kWh") + ylab("Frequency") +	ggtitle("Electricity") + xlim(0, 1050) + geom_text(data = annotations1, aes(x = x, y = y, label = paste(label, x)), size = 3, fontface = "bold")
fig2 <- ggplot(data = data, aes(WatPred)) + geom_histogram(bins = 40, color = "#000000", fill = "#0099F8") + 	xlab(TeX("$M^3$")) + ylab("Frequency") +	ggtitle("Water") + xlim(0, 100) + geom_text(data = annotations2, aes(x = x, y = y, label = paste(label, x)), size = 3, fontface = "bold")
fig3 <- ggplot(data = data, aes(GasPred)) + geom_histogram(bins = 40, color = "#000000", fill = "#0099F8") + 	xlab(TeX("$M^3$")) + ylab("Frequency") +	ggtitle("Gas") + xlim(0, 80) + geom_text(data = annotations3, aes(x = x, y = y, label = paste(label, x)), size = 3, fontface = "bold") 
# FIG <- ggarrange(fig1, fig2, fig3, ncol = 1, nrow = 3)
# FIG
# annotate_figure(FIG, top = text_grob("Utilities: Predictive distribution", color = "black", face = "bold", size = 14))

########################## Sleep in PhD students ########################## 
##### Beta as prior #########################
rm(list = ls())
set.seed(010101)
# Function to calibrate the hyperparameters of the beta prior distribution
elec <- function(theta, perc, q){
  E1 <- perc[1] - pbeta(q[1], theta[1], theta[2])
  E2 <- perc[2] - pbeta(q[2], theta[1], theta[2])
  loss <- E1^2 + E2^2
  return(loss)
}
q1 <- 0.4; q2 <- 0.75; q <- c(q1, q2) # Quantiles
p1e <- 0.6; p2e <- 0.95; perc <- c(p1e, p2e) # Cumulate probability from experts 
theta0 <- c(1, 1) # Parameters
elec(theta0, perc, q)
Res <- optim(theta0, elec, perc = perc, q = q)
Res$par
p <- seq(0, 1, length = 500)
a <- Res$par[1]
b <- Res$par[2]
pbeta(q1, a, b) # P(p <= q1) = p1e
pbeta(q2, a, b) # P(p <= q2) = p2e
s = 11
f = 16

prior=dbeta(p,a,b)
like=dbeta(p,s+1,f+1) #Seeing the binomial distribution as a likelihood is like a beta distribution!!!
post=dbeta(p,a+s,b+f)
plot(p,post,type="l",ylab="Density",lty=2,lwd=3,col=gray.colors(1,start=0.1))
lines(p,like,lty=1,lwd=3,col=gray.colors(1,start=0.4))
lines(p,prior,lty=3,lwd=3,col=gray.colors(1,start=0.6))
legend(.7,4,c("Prior","Likelihood","Posterior"),
       lty=c(3,1,2),lwd=c(3,3,3),col=c(col=gray.colors(1,start=0.6),col=gray.colors(1,start=0.4),col=gray.colors(1,start=0.1)))



p = c(seq(0.05, 0.95, by = 0.1))
prior = c(2, 4, 8, 8, 4, 2, 1, 1, 1, 1)
prior = prior/sum(prior)
data = c(11, 16)
####### Verosimilitu #################
s=11;f=16
#p1 = p + 0.5 * (p == 0) - 0.5 * (p == 1)
llike = s * log(p) + f * log(1 - p)
#like = p^s*(1-p)^f
#like1 = like * (p > 0) * (p < 1) - 999 * ((p == 0) * (s > 0) + (p == 1) * (f > 0))
like2 = exp(llike -max(llike)) #Escalar la likelihood para que se vea bien!!
product = like2 * prior
post = product/sum(product)

#pdisc(p,prior,data)

##############################################3
#post = pdisc(p, prior, data)
plot(p, prior, type = "l", ylab="Probability",lwd=2,lty=2,ylim=c(0,0.99),
     col=gray.colors(1,start=0.1))
lines(p,post,type = "l",lwd=2,col=gray.colors(1,start=0.6))
lines(p,like2,type = "l",lwd=4,col=gray.colors(1,start=0.3))
legend(0.7,0.85,c('Prior','Likelihood','Post'),lty=c(2,1,1),
       col=c(gray.colors(1,start=0.1),gray.colors(1,start=0.3),gray.colors(1,start=0.6)),lwd=2)

