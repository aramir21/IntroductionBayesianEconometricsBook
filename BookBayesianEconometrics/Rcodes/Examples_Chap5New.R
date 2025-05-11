########################## Mining change point ########################## 
# Clean workspace and set seed
rm(list = ls())
set.seed(10101)  # Avoid leading 0 (R interprets as octal)

# Load data
dataset <- read.csv(
  "https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/MiningDataCarlin.csv",
  header = TRUE
)

str(dataset)

# Hyperparameters
a10 <- 0.5; b10 <- 1
a20 <- 0.5; b20 <- 1

# Extract data
y <- dataset$Count
sum_y <- sum(y)
N <- length(y)

# Storage for posterior samples
S <- 10000
theta1 <- numeric(S)
theta2 <- numeric(S)
kk <- numeric(S)

# Initial value of change point
k <- 60

# Gibbs sampler
for (i in 1:S) {
  a1 <- a10 + sum(y[1:k])
  b1 <- b10 + k
  theta1[i] <- rgamma(1, shape = a1, rate = b1)
  
  a2 <- a20 + sum(y[(k + 1):N])
  b2 <- b20 + (N - k)
  theta2[i] <- rgamma(1, shape = a2, rate = b2)
  
  pp <- numeric(N)
  for (l in 1:N) {
    pp[l] <- exp(l * (theta2[i] - theta1[i])) *
      (theta1[i] / theta2[i])^sum(y[1:l])
  }
  
  prob <- pp / sum(pp)
  k <- sample(1:N, 1, prob = prob)
  kk[i] <- k
}

# Summaries
library(coda)
summary(mcmc(theta1))
summary(mcmc(theta2))
summary(mcmc(kk))

# Plot histogram
hist(
  kk,
  main = "Histogram: Posterior mean change point",
  xlab = "Posterior mean",
  col = "blue",
  breaks = 25
)
