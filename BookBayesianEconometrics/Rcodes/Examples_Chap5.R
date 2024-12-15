########################## Mining change point ########################## 
rm(list = ls())
set.seed(010101)
dataset<-read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/MiningDataCarlin.csv",header=T)
attach(dataset)
str(dataset)
a10 <- 0.5; a20 <- 0.5
b10 <- 1; b20 <- 1
y <- Count
sumy <- sum(Count); T <- length(Count)
theta1 <- NULL; theta2 <- NULL
kk <- NULL; H <- 60
MCMC <- 20000; burnin <- 1000; S <- MCMC + burnin; keep <- (burnin+1):S
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
for(s in 1:S){
  a1 <- a10 + sum(y[1:H])
  b1 <- b10+H
  theta11 <- rgamma(1,a1,b1)
  theta1 <- c(theta1,theta11)
  a2 <- a20 + sum(y[(1+H):T])
  b2 <- b20 + T-H
  theta22 <- rgamma(1,a2,b2)
  theta2 <- c(theta2,theta22)
  pp<-NULL
  for(l in 1:T){
    p <- exp(l*(theta22-theta11))*(theta11/theta22)^(sum(y[1:l]))
    pp <- c(pp,p)
  }
  prob <- pp/sum(pp)
  H <- sample(1:T,1,prob=prob)
  kk <- c(kk,H)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
library(coda)
theta1Post <- mcmc(theta1[keep])
plot(theta1Post); summary(theta1Post); autocorr.plot(theta1Post)
theta2Post <- mcmc(theta2[keep])
plot(theta2Post); summary(theta2Post); autocorr.plot(theta2Post)
HPost <- mcmc(kk); plot(HPost); summary(HPost); autocorr.plot(HPost)
hist(HPost, main = "Histogram: Posterior mean change point", xlab = "Posterior mean", col = "blue", breaks = 25)

########################## Metropolis-Hastings: Beta posterior########################## 
rm(list = ls())
set.seed(010101)
an <- 16.55; bn <- 39.57
S <- 100000
p <- runif(S)
acept <- rep(0, S)
for (s in 2:S){
  pc <- runif(1)
  a <- dbeta(pc, an, bn)/dbeta(p[s-1], an, bn)
  U <- runif(1)
  if(U <= a){
    p[s] <- pc
    acept[s] <- 1
  }else{
    p[s] <- p[s-1]
    acept[s] <- 0
  }
}
mean(acept)
mean(p); sd(p)
an/(an + bn); (an*bn/((an+bn)^2*(an+bn+1)))^0.5
h <- hist(p, breaks=50, col="blue", xlab="Proportion Ph.D. students sleeping at least 6 hours", main="Beta draws from a Metropolis-Hastings algorithm")
pfit <- seq(min(p),max(p),length=50)
yfit<-dbeta(pfit, an, bn)
yfit <- yfit*diff(h$mids[1:2])*length(p)
lines(pfit, yfit, col="red", lwd=2)

########################## Hamiltonian Monte Carlo: Mixture distribution ########################## 
rm(list = ls())
set.seed(010101)
N <- 5000
mu1 <- 1; mu2 <- -2
p1 <- 1/3; p2 <- 2/3
id <- sample(c(1,2), N, replace = TRUE, prob = c(p1, p2))
Y <- rep(NA, N)
for(i in 1:N){
  if(id[i] == 1){
    Y[i] <- rnorm(1, mu1)
  }else{
    Y[i] <- rnorm(1, mu2)
  }
}
plot(density(Y), main = "Mixture of normals", xlab = "Values")
loglikelihood <- function(mu1, mu2){
  sum(log((p1*dnorm(Y-mu1) + p2*dnorm(Y-mu2))))
}

# Create a grid of mu and sigma values
mu1s <- seq(-5, 5, length.out = 200)
mu2s <- seq(-5, 5, length.out = 200)

# Evaluate the log-likelihood over the grid
log_lik_matrix <- outer(mu1s, mu2s, 
                        Vectorize(function(mu1, mu2) loglikelihood(mu1, mu2)))

plot(x = 0, y = 0, type = "n", xlim = c(-5, 5), ylim = c(-5, 5), xlab = "", ylab = "")
u <- c(-5, 5, -5, 5)
tcol <- terrain.colors(12)
rect(u[1], u[3], u[2], u[4], col = tcol[8], border = "red")
contour(mu1s, mu2s, log_lik_matrix, 
        xlab = expression(mu1), ylab = expression(mu2), 
        main = "Log-Likelihood Contour Plot", add = TRUE)


da=rbind(rnorm(10^2),2.5+rnorm(3*10^2))
like=function(mu){
  sum(log((.25*dnorm(da-mu[1])+.75*dnorm(da-mu[2]))))}

scale=1
the=matrix(runif(2,-2,5),ncol=2)
curlike=hval=like(x)
Niter=10^4
for (iter in (1:Niter)){
  prop=the[iter,]+rnorm(2)*scale
  if ((max(-prop)>2)||(max(prop)>5)||
      (log(runif(1))>like(prop)-curlike)) prop=the[iter,]
  curlike=like(prop)
  hval=c(hval,curlike)
  the=rbind(the,prop)}
plot(density(the))
dim(the)



# Load necessary libraries
library(MASS)
library(ggplot2)

# Step 1: Simulate data from a mixture of two normal distributions
set.seed(42)
n <- 500
mu1 <- 1    # Mean of the first population
mu2 <- 3    # Mean of the second population
sigma <- 1  # Common standard deviation
p <- 0.5    # Proportion of each sub-population

# Generate data
z <- rbinom(n, 1, p)  # Latent variable indicating which component
data <- z * rnorm(n, mean = mu1, sd = sigma) +
  (1 - z) * rnorm(n, mean = mu2, sd = sigma)

plot(density(data))
abline(v=mu1); abline(v=mu2); abline(v=(mu1+mu2)/2, col = "red")

# Step 2: Define the Bayesian model
prior_mean <- 0
prior_sd <- 10

# Log-posterior function to avoid underflow
log_posterior <- function(mean) {
  log_likelihood <- sum(dnorm(data, mean = mean, sd = sigma, log = TRUE))
  log_prior <- dnorm(mean, mean = prior_mean, sd = prior_sd, log = TRUE)
  return(log_likelihood + log_prior)
}

# Step 3: Evaluate the posterior distribution
grid <- seq(-1, 11, length.out = 1000)  # Grid of potential mean values
log_posterior_values <- sapply(grid, log_posterior)

# Convert log-posterior to posterior density
posterior_values <- exp(log_posterior_values - max(log_posterior_values))  # Scale to avoid overflow
posterior_values <- posterior_values / sum(posterior_values * diff(grid)[1])  # Normalize

# Step 4: Plot the posterior distribution
posterior_df <- data.frame(mean = grid, density = posterior_values)

ggplot(posterior_df, aes(x = mean, y = density)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Posterior Distribution of the Mean",
       x = "Mean",
       y = "Density") +
  theme_minimal()

# Step 5: Overlay the histogram of data for visualization
ggplot() +
  geom_histogram(aes(x = data, y = ..density..), bins = 30, fill = "gray", alpha = 0.6) +
  geom_line(data = posterior_df, aes(x = mean, y = density), color = "blue", size = 1) +
  labs(title = "Histogram of Data and Posterior Distribution",
       x = "Value",
       y = "Density") +
  theme_minimal()

