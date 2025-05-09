########### Chapter 4 ###############

########### Dirichlet-Multinomial model: Liverpool vs Manchester City ##############

# Clear workspace and set reproducible seed
rm(list = ls())
set.seed(10101)

########################## Multinomial-Dirichlet example: Liverpool vs Manchester City ##########################

# Load required libraries
library(dplyr)
library(readr)
library(sirt)
library(MCMCpack)
library(Compositional)

# Load data
data <- read_csv(
  "https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/DataOddsLIVvsMAN.csv",
  show_col_types = FALSE
)

# Match and odds info:
# - Odds collected: 2022-10-06 19:00 (Colombia time)
# - Match played: 2022-10-16
# Source: https://www.oddsportal.com/soccer/england/premier-league/liverpool-manchester-city-WrqgEz5S/

# Compute implied probabilities from odds
probs <- data %>%
  mutate(
    pns1 = 1 / home,
    pns2 = 1 / draw,
    pns3 = 1 / away,
    sum_inv_odds = pns1 + pns2 + pns3,
    p1 = pns1 / sum_inv_odds,
    p2 = pns2 / sum_inv_odds,
    p3 = pns3 / sum_inv_odds
  ) %>%
  dplyr::select(p1, p2, p3)

# Estimate Dirichlet parameters from betting odds
dir_mle <- dirichlet.mle(probs)
alpha0_odds <- dir_mle$alpha
print(alpha0_odds)

# Historical results: 2 wins Liverpool, 2 draws, 1 win Man City
# Source: https://www.11v11.com/teams/manchester-city/tab/opposingTeams/opposition/Liverpool/
y <- c(2, 2, 1)

# Marginal log-likelihood function (negative for optimization)
marginal_likelihood <- function(alpha) {
  n <- sum(y)
  res1 <- sum(sapply(seq_along(y), function(i) {
    lgamma(alpha[i] + y[i]) - lgamma(alpha[i])
  }))
  res <- lgamma(sum(alpha)) - lgamma(sum(alpha) + n) + res1
  return(-res)
}

# Empirical Bayes estimate using MLE
emp_bayes <- optim(alpha0_odds, marginal_likelihood, method = "BFGS", control = list(maxit = 10000))
alpha0_eb <- emp_bayes$par
print(alpha0_eb)

# Compare with prior based on odds
print(alpha0_odds)

# Bayes factor: empirical Bayes vs odds-based prior
bf <- exp(-marginal_likelihood(alpha0_eb)) / exp(-marginal_likelihood(alpha0_odds))
print(bf)

# Posterior parameters (empirical Bayes prior + data)
alpha_n <- alpha0_eb + y

# Draw posterior samples
set.seed(10101)
S <- 100000
thetas <- rdirichlet(S, alpha_n)
colnames(thetas) <- c("Liverpool", "Draw", "Manchester")
head(thetas)

# Visualize posterior distribution
bivt.contour(thetas, cont.line = FALSE, appear = FALSE)

# Predictive distribution for hypothetical outcome
y0 <- c(2, 0, 3)
pred_draws <- apply(thetas, 1, function(p) {
  rmultinom(1, size = sum(y0), prob = p)
})

prob_y0_sim <- mean(apply(pred_draws, 2, function(draw) all(draw == y0)))
print(prob_y0_sim)

# Predictive probability using closed-form expression
predictive_prob_y0 <- function(y0) {
  n <- sum(y0)
  res1 <- sum(sapply(seq_along(y), function(i) {
    lgamma(alpha_n[i] + y0[i]) - lgamma(alpha_n[i]) - lfactorial(y0[i])
  }))
  res <- lfactorial(n) + lgamma(sum(alpha_n)) - lgamma(sum(alpha_n) + n) + res1
  return(exp(res))
}

predictive_prob_y0(y0)

########################## The multivariate normal-normal/inverse-Wishart model example ########################## 
# Tangent portfolio

# Load required libraries
library(quantmod)
library(xts)
library(ggplot2)
library(gridExtra)
library(purrr)
library(dplyr)

# Define date range
start_date <- as.Date("2021-01-01")
end_date <- as.Date("2022-09-30")
dates <- seq(start_date, end_date, by = "day")

# Tickers of interest
tickers <- c("AAPL", "NFLX", "AMZN", "GOOG", "INTC", "META", "MSFT", "TSLA", "NVDA", "PYPL")
p <- length(tickers)

# Download adjusted closing prices
getSymbols(tickers, from = start_date, to = end_date, auto.assign = TRUE)

prices <- map(tickers, ~ Ad(get(.x))) %>%
  reduce(merge) %>%
  `colnames<-`(tickers) %>%
  as.data.frame()

# Calculate daily log returns
returns <- apply(prices, 2, function(x) diff(log(x))) %>%
  as.data.frame()

# Download 10-year Treasury yield from FRED
t10yr <- getSymbols("DGS10", src = "FRED", from = start_date, to = end_date, auto.assign = FALSE)
t10yr_daily <- ((1 + t10yr / 100)^(1 / 365)) - 1
t10yr_daily <- t10yr_daily[rownames(returns), ]

# Compute excess returns
excess_returns <- as.matrix(returns) - kronecker(t(rep(1, p)), as.matrix(t10yr_daily))

# Convert to data frame with dates
df <- as.data.frame(excess_returns)
df$Date <- as.Date(rownames(df))
df$Month <- months(df$Date)
df$Year <- format(df$Date, "%y")

# Aggregate monthly means
monthly_means <- map(1:p, function(i) {
  aggregate(df[[i]] ~ Month + Year, data = df, FUN = mean)
})

# Extract values into matrix
data_excess <- matrix(0, nrow = nrow(monthly_means[[1]]), ncol = p)
for (i in 1:p) {
  data_excess[, i] <- as.numeric(monthly_means[[i]][, 3])
}
colnames(data_excess) <- tickers

# Hyperparameters
N <- nrow(data_excess)
mu_0 <- rep(0, p)
beta_0 <- 1
psi_0 <- 100 * diag(p)
alpha_0 <- p + 2

# Posterior parameters
alpha_n <- N + alpha_0
v_n <- alpha_n + 1 - p
mu_hat <- colMeans(data_excess)
mu_n <- (N / (N + beta_0)) * mu_hat + (beta_0 / (N + beta_0)) * mu_0
S <- t(data_excess - matrix(mu_hat, N, p, byrow = TRUE)) %*% 
  (data_excess - matrix(mu_hat, N, p, byrow = TRUE))
psi_n <- psi_0 + S + (N * beta_0 / (N + beta_0)) * 
  tcrossprod(mu_hat - mu_0)

beta_n <- N + beta_0
sigma_n <- psi_n / ((alpha_n + 1 - p) * beta_n)
cov_n <- (sigma_n * (1 + beta_n)) * v_n / (v_n - 2)
cov_inv <- solve(cov_n)

# Optimal portfolio weights (Bayesian mean-variance)
opt_weights <- t(cov_inv %*% mu_n / as.numeric(t(rep(1, p)) %*% cov_inv %*% mu_n))
colnames(opt_weights) <- tickers

# Result
opt_weights
