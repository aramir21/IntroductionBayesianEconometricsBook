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
