########################## Math test example ########################## 
rm(list = ls())
set.seed(010101)
N <- 50
# Sample size
muhat <- 102
# Sample mean
sig2hat <- 100
# Sample variance

# Hyperparameters
mu0 <- 100
sig20 <- 100
delta0 <- 0.001
alpha0 <- 0.001

MCMC <- 10000; burnin <- 1000; S <- MCMC + burnin
keep <- (burnin+1):S
# Posterior draws
alphan <- alpha0 + N
sig2Post <- rep(NA, S)
muPost <- rep(NA, S)
sig2 <- sig20
for(s in 1:S){
  sig2n <- (1/(sig2/N)+1/sig20)^(-1)
  mun <- sig2n*(muhat/(sig2/N)+mu0/sig20)
  mu <- rnorm(1, mun, sig2n^0.5)
  muPost[s] <- mu
  deltan <- N*(sig2hat + (muhat - mu)^2)
  sig2 <- invgamma::rinvgamma(1, shape = alphan, rate = deltan)
  sig2Post[s] <- sig2
}
sig2s <- coda::mcmc(sig2Post[keep]) 
mus <- coda::mcmc(muPost[keep]) 
summary(sig2s)
summary(mus)
hist(mus, main = "Histogram: Posterior mean", xlab = "Posterior mean", col = "blue", breaks = 50)
muPost_tq <- quantile(mus, c(0.025, 0.5, 0.975))
muPost_tq
cutoff <- 103
PmuPost_tcutoff <- mean(mus > cutoff)
PmuPost_tcutoff

########################## Metropolis-Hastings: Cauchy distribution ########################## 
rm(list = ls())
set.seed(010101)
S <- 100000
df <- 5
theta1 <- runif(S)
acept1 <- rep(0, S)
theta2 <- runif(S)
acept2 <- rep(0, S)
for (s in 2:S){
  thetac1 <- rnorm(1)
  thetac2 <- rt(1, df = df)
  a1 <- (dcauchy(thetac1)*dnorm(theta1[s-1]))/(dcauchy(theta1[s-1])*dnorm(thetac1))
  a2 <- (dcauchy(thetac2)*dt(theta2[s-1], df = df))/(dcauchy(theta2[s-1])*dt(thetac2, df = df))
  U1 <- runif(1)
  if(U1 <= a1){
    theta1[s] <- thetac1
    acept1[s] <- 1
  }else{
    theta1[s] <- theta1[s-1]
    acept1[s] <- 0
  }
  if(U1 <= a2){
    theta2[s] <- thetac2
    acept2[s] <- 1
  }else{

    theta2[s] <- theta2[s-1]
    acept2[s] <- 0
  }
}
mean(acept1); mean(acept2)
mean(theta1); sd(theta1)
mean(theta2); sd(theta2)
plot(coda::mcmc(theta1)); coda::autocorr.plot(coda::mcmc(theta1))
plot(coda::mcmc(theta2)); coda::autocorr.plot(coda::mcmc(theta2))
h <- hist(theta1, breaks=50, col="blue", xlab="x", main="Cauchy draws from a Metropolis-Hastings algorithm: Normal standard proposal")
pfit <- seq(min(theta1),max(theta1),length=50)
yfit<-dcauchy(pfit)
yfit <- yfit*diff(h$mids[1:2])*length(theta1)
lines(pfit, yfit, col="red", lwd=2)

h <- hist(theta2, breaks=50, col="blue", xlab="x", main="Cauchy draws from a Metropolis-Hastings algorithm: Student's t proposal")
pfit <- seq(min(theta2),max(theta2),length=50)
yfit<-dcauchy(pfit)
yfit <- yfit*diff(h$mids[1:2])*length(theta2)
lines(pfit, yfit, col="red", lwd=2)
