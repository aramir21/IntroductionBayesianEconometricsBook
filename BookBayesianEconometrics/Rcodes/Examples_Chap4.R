rm(list = ls())
set.seed(010101)
########################## Multinomial-Dirichlet example: Liverpool vs Manchester city ##########################
Data <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/DataOddsLIVvsMAN.csv", sep = ",", header = TRUE, quote = "")
# https://www.oddsportal.com/soccer/england/premier-league/liverpool-manchester-city-WrqgEz5S/
# Betting odds at 6-10-2022 19:00 hours (Colombia time)
# Liverpool vs Manchester city at 16-10-2022
attach(Data)
library(dplyr)
Probs <- Data %>%
  mutate(pns1 = 1/home, pns2 = 1/draw, pns3 = 1/away) %>%
  mutate(SumInvOdds = pns1 + pns2 + pns3) %>% 
  mutate(p1 = pns1/SumInvOdds, p2 = pns2/SumInvOdds, p3 = pns3/SumInvOdds) %>% 
  select(p1, p2, p3)

DirMLE <- sirt::dirichlet.mle(Probs)
alpha0odds <- DirMLE$alpha
alpha0odds
y <- c(2, 2, 1) 
# Liverpool win, draw, Manchester city win
# https://www.11v11.com/teams/manchester-city/tab/opposingTeams/opposition/Liverpool/
# Last 5 matches in Premier League Liverpool local: 14/01/2018-10/04/2022

# Marginal likelihood
MarLik <- function(a0){
  n <- sum(y)
  Res1 <- sum(sapply(1:length(y), function(l){lgamma(a0[l] + y[l]) - lgamma(a0[l])}))
  Res <- lgamma(sum(a0)) - lgamma(sum(a0)+n) + Res1
  return(-Res)
}

EmpBay <- optim(alpha0odds, MarLik, method = "BFGS", control = list(maxit = 10000))
alpha0EB <- EmpBay$par
alpha0EB
alpha0odds

BF <- exp(-MarLik(alpha0EB))/exp(-MarLik(alpha0odds))
BF

set.seed(010101)
alphan <- alpha0odds + y 
S <- 100000 
# Sample draws
thetas <- MCMCpack::rdirichlet(S, alphan)
colnames(thetas) <- c("Liverpool", "Draw", "Manchester")
head(thetas)
Compositional::bivt.contour(thetas, cont.line = FALSE, appear = FALSE)

# Predictive
y0 <- c(2, 0, 3)
Pred <- apply(thetas, 1, function(p) {rmultinom(1, size = sum(y0), prob = p)})
ProbY0 <- sum(sapply(1:S, function(s) {sum(Pred[,s] == y0) == 3}))/S
ProbY0
PredY0 <- function(y0){
  n <- sum(y0)
  Res1 <- sum(sapply(1:length(y), function(l){lgamma(alphan[l] + y0[l]) - lgamma(alphan[l]) - lfactorial(y0[l])}))
  Res <- lfactorial(n) + lgamma(sum(alphan)) - lgamma(sum(alphan)+n) + Res1
  return(exp(Res))
}
PredY0(y0) 

########################## The multivariate normal-normal/inverse-Wishart model example ########################## 
# Tangency portfolio

library(quantmod)
library(xts)
library(ggplot2)
library(gridExtra) # grid.arrange

graphics.off()
rm(list=ls())

# Data Range
sdate <- as.Date("2021-01-01")
edate <- as.Date("2022-09-30")
Date <- seq(sdate, edate, by = "day")

tickers <- c("AAPL", "NFLX", "AMZN", "GOOG", "INTC",
             "META", "MSFT", "TSLA", "NVDA", "PYPL")
p <- length(tickers)

# AAPL: Apple, NFLX: Netflix, AMZN: Amazon, MSFT: Microsoft
# GOOG: Google, META: Meta, TSLA: Tesla, NVDA: NVIDIA Corporation
# INTC: Intel, PYPL: PayPal 
ss_stock <- getSymbols(tickers, from=sdate, to=edate, auto.assign = T)
ss_stock <- purrr::map(tickers,function(x) Ad(get(x)))
ss_stock <- as.data.frame(purrr::reduce(ss_stock, merge))
colnames(ss_stock) <- tickers
ss_rtn <- as.data.frame(apply(ss_stock, 2, function(x) {diff(log(x), 1)}))
# Daily returns

t10yr <- getSymbols(Symbols = "DGS10", src = "FRED", from=sdate, to=edate, auto.assign = F)
# 10-Year US Treasury yield data from the Federal Reserve Electronic Database (FRED)
t10yrd <- (1 + t10yr/100)^(1/365)-1 
# Daily returns
t10yrd <- t10yrd[row.names(ss_rtn)]
Exc_rtn <- as.matrix(ss_rtn) - kronecker(t(rep(1, p)), as.matrix(t10yrd))
# Excesses of return
df <- as.data.frame(Exc_rtn)
df$Date <- as.Date(rownames(df))
#  Get months
df$Month <- months(df$Date)
#  Get years
df$Year <- format(df$Date, format="%y")
head(df)

#  Aggregate on months and year and get mean
Data <- sapply(1:p, function(i) {aggregate(df[, i] ~ Month + Year, df, mean)})
DataExcRtn <- matrix(0, length(Data[, 1]$Month), p)
for(i in 1:p){
  DataExcRtn[, i] <- as.numeric(Data[, i]$`df[, i]`)
}

colnames(DataExcRtn) <- tickers
head(DataExcRtn)
# draw graph
plot1 <- ggplot(df, aes(x = Date, y = AAPL)) +
  geom_line(color = "blue", linewidth=1.2) + 
  ggtitle("Apple excess of return") + xlab("Date") + ylab("Excess return") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_date(date_labels = "%y-%m", date_breaks = "2 months")

plot2 <- ggplot(df, aes(x = Date, y = NFLX)) +
  geom_line(color = "blue", linewidth=1.2) + 
  ggtitle("Netflix excess of return") + xlab("Date") + ylab("Excess return") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_date(date_labels = "%y-%m", date_breaks = "2 months")

grid.arrange(plot1, plot2, ncol=2, nrow = 1)

# Hyperparameters #
N <- dim(DataExcRtn)[1]
mu0 <- rep(0, p)
beta0 <- 1
Psi0 <- 100 * diag(p)
alpha0 <- p + 2

# Posterior parameters #
alphan <- N + alpha0
vn <- alphan + 1 - p
muhat <- colMeans(DataExcRtn)
mun <- N/(N + beta0) * muhat + beta0/(N + beta0) * mu0
S <- t(DataExcRtn - rep(1, N)%*%t(muhat))%*%(DataExcRtn - rep(1, N)%*%t(muhat)) 
Psin <- Psi0 + S + N*beta0/(N + beta0)*(muhat - mu0)%*%t(muhat - mu0)
betan <- N + beta0
Sigman <- Psin/((alphan + 1 - p)*betan)
Covarn <- (Sigman * (1 + betan)) * vn / (vn - 2)
Covari <- solve(Covarn)
OptShare <- t(Covari%*%mun/as.numeric((t(rep(1, p))%*%Covari%*%mun)))
colnames(OptShare) <- tickers
OptShare

########################## Linear regression example ##########################
rm(list = ls())
set.seed(010101)
# Electricity demand
DataUt <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/Utilities.csv", sep = ",", header = TRUE, quote = "")

DataUtEst <- DataUt %>% 
  filter(Electricity	 != 0)

attach(DataUtEst)
# Dependent variable: Monthly consumption (kWh) in log
Y <- log(Electricity) 

# Regressors quantity including intercept
X <- cbind(LnPriceElect, IndSocio1, IndSocio2, Altitude, Nrooms, HouseholdMem, Children, 
          Lnincome, 1)
# LnPriceElect: Price per kWh (USD) in log
# IndSocio1, IndSocio2, IndSocio3: Indicators socio-economic condition (1) is the lowest and (3) the highest
# Altitude: Indicator of household location (1 is more than 1000 meters above sea level)
# Nrooms: Number of rooms in house
# HouseholdMem: Number of household members
# Children: Indicator por presence of children in household (1)
# Lnincome: Monthly income (USD) in log
k <- dim(X)[2]
N <- dim(X)[1]

# Hyperparameters
d0 <- 0.001/2
a0 <- 0.001/2
b0 <- rep(0, k)
c0 <- 1000
B0 <- c0*diag(k)

# Posterior parameters
bhat <- solve(t(X)%*%X)%*%t(X)%*%Y
Bn <- as.matrix(Matrix::forceSymmetric(solve(solve(B0) + t(X)%*%X))) # Force this matrix to be symmetric
bn <- Bn%*%(solve(B0)%*%b0 + t(X)%*%X%*%bhat)
dn <- as.numeric(d0 + t(Y)%*%Y+t(b0)%*%solve(B0)%*%b0-t(bn)%*%solve(Bn)%*%bn)
an <- a0 + N
Hn <- Bn*dn/an


# Posterior draws
S <- 10000 # Number of draws from posterior distributions
sig2 <- MCMCpack::rinvgamma(S,an/2,dn/2)
summary(coda::mcmc(sig2))
Betas <- LaplacesDemon::rmvt(S, bn, Hn, an)
summary(coda::mcmc(Betas))

# Log marginal function (multiply by -1 due to minimization)
LogMarLikLM <- function(X, c0){
  k <- dim(X)[2]
  N <- dim(X)[1]
  
  # Hyperparameters
  B0 <- c0*diag(k)
  b0 <- rep(0, k)
  
  # Posterior parameters
  bhat <- solve(t(X)%*%X)%*%t(X)%*%Y
  Bn <- as.matrix(Matrix::forceSymmetric(solve(solve(B0) + t(X)%*%X))) # Force this matrix to be symmetric
  bn <- Bn%*%(solve(B0)%*%b0 + t(X)%*%X%*%bhat)
  dn <- as.numeric(d0 + t(Y)%*%Y+t(b0)%*%solve(B0)%*%b0-t(bn)%*%solve(Bn)%*%bn)
  an <- a0 + N
  # Log marginal likelihood
  logpy <- (N/2)*log(1/pi)+(a0/2)*log(d0)-(an/2)*log(dn) + 0.5*log(det(Bn)/det(B0)) + lgamma(an/2)-lgamma(a0/2)
  
  return(-logpy)
}

cs <- c(10^0, 10^3, 10^6, 10^10, 10^12, 10^15, 10^20)
# Observe -1 to recover the right sign
LogML <- sapply(cs, function(c) {-LogMarLikLM(c0=c, X = X)}) 

# Regressor without price
Xnew <- cbind(IndSocio1, IndSocio2, Altitude, Nrooms, HouseholdMem, Children, 
              Lnincome, 1)
# Observe -1 to recover the right sign
LogMLnew <- sapply(cs, function(c) {-LogMarLikLM(c0=c, X = Xnew)})

# Bayes factor
BF <- exp(LogML - LogMLnew)
BF

# Predictive distribution
Xpred <- c(log(0.15), 1, 0, 0, 2, 3, 1, log(500), 1)
Mean <- Xpred%*%bn
Hn <- dn*(1+t(Xpred)%*%Bn%*%Xpred)/an
ExpKwH <- exp(LaplacesDemon::rmvt(S, Mean, Hn, an))
summary(ExpKwH)
HDI <- HDInterval::hdi(ExpKwH, credMass = 0.95) # Highest posterior density credible interval
HDI
hist(ExpKwH, main = "Histogram: Monthly demand of electricity", xlab = "Monthly kWh", col = "blue", breaks = 50)
