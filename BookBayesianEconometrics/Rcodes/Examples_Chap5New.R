########################## Mining change point ########################## 
# Clean workspace and set seed
rm(list = ls())
set.seed(10101)  

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

########################## Metropolis-Hastings: Beta posterior########################## 
# Clear workspace and set seed for reproducibility
rm(list = ls())
set.seed(10101)

# Beta distribution parameters
a_n <- 16.55
b_n <- 39.57

# Number of samples
S <- 100000

# Initialize vectors
samples <- numeric(S)
accept_flags <- logical(S)

# Initial value
samples[1] <- runif(1)

# Metropolis-Hastings sampling
for (s in 2:S) {
  candidate <- runif(1)
  
  alpha <- dbeta(candidate, a_n, b_n) / dbeta(samples[s - 1], a_n, b_n)
  u <- runif(1)
  
  if (u <= alpha) {
    samples[s] <- candidate
    accept_flags[s] <- TRUE
  } else {
    samples[s] <- samples[s - 1]
    accept_flags[s] <- FALSE
  }
}

# Posterior summaries
mean_accept <- mean(accept_flags)
mean_sample <- mean(samples)
sd_sample <- sd(samples)

# True mean and standard deviation of Beta(a_n, b_n)
true_mean <- a_n / (a_n + b_n)
true_sd <- sqrt((a_n * b_n) / ((a_n + b_n)^2 * (a_n + b_n + 1)))

# Print summaries
cat("Acceptance rate:", mean_accept, "\n")
cat("Sample mean:", mean_sample, " | True mean:", true_mean, "\n")
cat("Sample SD:  ", sd_sample, " | True SD:  ", true_sd, "\n")

# Histogram and overlaid density
hist_info <- hist(
  samples,
  breaks = 50,
  col = "blue",
  xlab = "Proportion of Ph.D. students sleeping ≥ 6 hours",
  main = "Beta draws from a Metropolis-Hastings algorithm"
)

x_vals <- seq(min(samples), max(samples), length.out = 50)
y_vals <- dbeta(x_vals, a_n, b_n)
y_vals <- y_vals * diff(hist_info$mids[1:2]) * length(samples)

lines(x_vals, y_vals, col = "red", lwd = 2)

########################## Hamiltonian Monte Carlo: Bi-variate normal distribution ########################## 

# Load packages
library(MASS)
library(mvtnorm)
library(coda)
library(ggplot2)
library(ggpubr)
library(latex2exp)

set.seed(10101)

###----------------------------------------------
### Functions
###----------------------------------------------

# Gibbs sampler
gibbs_sampler <- function(theta, rho) {
  rnorm(1, mean = rho * theta, sd = sqrt(1 - rho^2))
}

# Metropolis-Hastings sampler
mh_sampler <- function(theta, rho, sig2) {
  sigma_target <- matrix(c(1, rho, rho, 1), 2, 2)
  sigma_proposal <- matrix(c(1, sig2, sig2, 1), 2, 2)
  
  theta_candidate <- mvrnorm(1, mu = theta, Sigma = sigma_proposal)
  
  a <- dmvnorm(theta_candidate, mean = c(0, 0), sigma = sigma_target) /
    dmvnorm(theta, mean = c(0, 0), sigma = sigma_target)
  
  if (runif(1) <= a) {
    list(theta = theta_candidate, accept = 1)
  } else {
    list(theta = theta, accept = 0)
  }
}

# Hamiltonian Monte Carlo
hmc_sampler <- function(theta, rho, epsilon, M) {
  sigma_target <- matrix(c(1, rho, rho, 1), 2, 2)
  L <- ceiling(1 / epsilon)
  M_inv <- solve(M)
  K <- length(theta)
  
  # Initial momentum
  momentum <- t(rmvnorm(1, mean = rep(0, K), sigma = M))
  theta_old <- theta
  
  # Initial log posterior + kinetic energy
  log_start <- dmvnorm(theta_old, mean = rep(0, K), sigma = sigma_target, log = TRUE) +
    dmvnorm(as.vector(momentum), mean = rep(0, K), sigma = M, log = TRUE)
  
  # Leapfrog steps
  for (l in 1:L) {
    grad <- -solve(sigma_target) %*% theta
    if (l == 1 || l == L) {
      momentum <- momentum + 0.5 * epsilon * grad
    } else {
      momentum <- momentum + epsilon * grad
    }
    theta <- theta + epsilon * M_inv %*% momentum
  }
  
  # Final log posterior + kinetic energy
  log_end <- dmvnorm(as.vector(theta), mean = rep(0, K), sigma = sigma_target, log = TRUE) +
    dmvnorm(as.vector(momentum), mean = rep(0, K), sigma = M, log = TRUE)
  
  alpha <- min(1, exp(log_end - log_start))
  
  if (runif(1) <= alpha) {
    theta_new <- as.vector(theta)
  } else {
    theta_new <- theta_old
  }
  
  list(theta = theta_new, prob = alpha)
}

###----------------------------------------------
### Parameters and storage
###----------------------------------------------

rho <- 0.98
sig2 <- 0.18^2
S_gibbs <- 8000
thin <- 20
K <- 2
keep <- seq(0, S_gibbs, by = thin)[-1]

theta_post_gibbs <- matrix(NA, S_gibbs, K)
theta_post_mh <- matrix(NA, S_gibbs, K)
accept_mh <- logical(S_gibbs)

theta_gibbs <- c(-2, 3)
theta_mh <- c(-2, 3)

###----------------------------------------------
### Gibbs and Metropolis-Hastings
###----------------------------------------------

for (s in 1:S_gibbs) {
  theta1 <- gibbs_sampler(theta_gibbs[2], rho)
  theta2 <- gibbs_sampler(theta1, rho)
  theta_gibbs <- c(theta1, theta2)
  
  res_mh <- mh_sampler(theta_mh, rho, sig2)
  theta_mh <- res_mh$theta
  
  theta_post_gibbs[s, ] <- theta_gibbs
  theta_post_mh[s, ] <- theta_mh
  accept_mh[s] <- res_mh$accept
}

cat("MH acceptance rate:", mean(accept_mh[keep]), "\n")

mcmc_gibbs <- mcmc(theta_post_gibbs[keep, ])
mcmc_mh <- mcmc(theta_post_mh[keep, ])

summary(mcmc_gibbs)
autocorr.plot(mcmc_gibbs)
plot(mcmc_mh)
autocorr.plot(mcmc_mh)

###----------------------------------------------
### Hamiltonian Monte Carlo
###----------------------------------------------

S_hmc <- 400
epsilon <- 0.05
M <- diag(2)

theta_post_hmc <- matrix(NA, S_hmc, K)
accept_prob_hmc <- numeric(S_hmc)
theta_hmc <- c(-2, 3)

for (s in 1:S_hmc) {
  res_hmc <- hmc_sampler(theta_hmc, rho, epsilon, M)
  theta_hmc <- res_hmc$theta
  theta_post_hmc[s, ] <- theta_hmc
  accept_prob_hmc[s] <- res_hmc$prob
}

mcmc_hmc <- mcmc(theta_post_hmc)
summary(mcmc_hmc)
summary(accept_prob_hmc)
plot(mcmc_hmc)
autocorr.plot(mcmc_hmc)

###----------------------------------------------
### Comparison Plot
###----------------------------------------------

df_plot <- data.frame(
  Iter = 1:S_hmc,
  HMC = theta_post_hmc[, 1],
  MH = theta_post_mh[keep, 1],
  Gibbs = theta_post_gibbs[keep, 1]
)

g1 <- ggplot(df_plot, aes(x = Iter)) +
  geom_point(aes(y = HMC), color = "black") +
  labs(x = "Iteration", y = TeX("$\\theta_{1}$"), title = "HMC algorithm")

g2 <- ggplot(df_plot, aes(x = Iter)) +
  geom_point(aes(y = MH), color = "black") +
  labs(x = "Iteration", y = TeX("$\\theta_{1}$"), title = "M-H algorithm")

g3 <- ggplot(df_plot, aes(x = Iter)) +
  geom_point(aes(y = Gibbs), color = "black") +
  labs(x = "Iteration", y = TeX("$\\theta_{1}$"), title = "Gibbs sampling")

ggarrange(g3, g2, g1, labels = c("A", "B", "C"), ncol = 3)

########################## Importance sampling: Cauchy distribution ########################## 
# Load required libraries
library(ggplot2)
library(ggpubr)
library(dplyr)

set.seed(10101)

# Parameters
S <- 20000  # Proposal sample size
L <- 10000  # Posterior sample size

###----------------------------------
### Importance Sampling: Normal proposal
###----------------------------------
theta_norm <- rnorm(S)
weights_norm <- dcauchy(theta_norm) / dnorm(theta_norm)
weights_norm <- weights_norm / sum(weights_norm)

theta_cauchy_norm <- sample(theta_norm, size = L, replace = TRUE, prob = weights_norm)
df_norm <- data.frame(x = theta_cauchy_norm)

g1 <- ggplot(df_norm, aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "skyblue", color = "black") +
  stat_function(fun = dcauchy, color = "red", linewidth = 1.2) +
  labs(title = "Cauchy draws via Importance Sampling",
       subtitle = "Normal(0,1) proposal",
       x = "x", y = "Density") +
  theme_minimal()

###----------------------------------
### Importance Sampling: Student's t proposal
###----------------------------------
df_t <- 3
theta_t <- rt(S, df = df_t)
weights_t <- dcauchy(theta_t) / dt(theta_t, df = df_t)
weights_t <- weights_t / sum(weights_t)

theta_cauchy_t <- sample(theta_t, size = L, replace = TRUE, prob = weights_t)
df_tpost <- data.frame(x = theta_cauchy_t)

g2 <- ggplot(df_tpost, aes(x = x)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "lightgreen", color = "black") +
  stat_function(fun = dcauchy, color = "red", linewidth = 1.2) +
  labs(title = "Cauchy draws via Importance Sampling",
       subtitle = "Student's t (df = 3) proposal",
       x = "x", y = "Density") +
  theme_minimal()

###----------------------------------
### Plot Importance Weights
###----------------------------------
df_weights <- data.frame(
  index = 1:S,
  Normal = weights_norm,
  Student_t = weights_t
)

g3 <- ggplot(df_weights) +
  geom_point(aes(x = index, y = Normal), color = "black", alpha = 0.6, size = 0.7) +
  geom_point(aes(x = index, y = Student_t), color = "blue", alpha = 0.6, size = 0.7) +
  labs(title = "Importance Sampling Weights",
       x = "Index", y = "Weight") +
  scale_y_continuous(trans = "log1p") +
  theme_minimal()

###----------------------------------
### Arrange and display plots
###----------------------------------
ggarrange(g1, g2, g3, ncol = 1, labels = c("A", "B", "C"))

########################## Sequential importance sampling: Dynamic linear model ########################## 

# Load packages
library(dplyr)
library(ggplot2)
library(latex2exp)

set.seed(10101)

# Parameters
n_particles <- 50000
sigma_w <- 1       # State noise
sigma_mu <- 1      # Observation noise
phi <- 0.5         # Observation coefficient
T <- 50            # Time points

###-------------------------------------
### Simulate true states and observations
###-------------------------------------
theta_true <- numeric(T)
theta_true[1] <- rnorm(1, mean = 0, sd = sigma_w)
for (t in 2:T) {
  theta_true[t] <- rnorm(1, mean = theta_true[t - 1], sd = sigma_w)
}
y_obs <- rnorm(T, mean = phi * theta_true, sd = sigma_mu)

###-------------------------------------
### Sequential Importance Sampling (SIS)
###-------------------------------------
particles <- matrix(0, nrow = n_particles, ncol = T)
weights <- matrix(0, nrow = n_particles, ncol = T)
weights_normalized <- matrix(0, nrow = n_particles, ncol = T)

# Initialization
particles[, 1] <- rnorm(n_particles, mean = 0, sd = sigma_w)
weights[, 1] <- dnorm(y_obs[1], mean = phi * particles[, 1], sd = sigma_mu)
weights_normalized[, 1] <- weights[, 1] / sum(weights[, 1])

# Recursive updates
for (t in 2:T) {
  particles[, t] <- rnorm(n_particles, mean = particles[, t - 1], sd = sigma_w)
  weights[, t] <- weights_normalized[, t - 1] * dnorm(y_obs[t], mean = phi * particles[, t], sd = sigma_mu)
  weights_normalized[, t] <- weights[, t] / sum(weights[, t])
}

###-------------------------------------
### Posterior filtering mean and std dev
###-------------------------------------
posterior_mean <- colSums(particles * weights_normalized)
posterior_sd <- sqrt(colSums((particles^2) * weights_normalized) - posterior_mean^2)

###-------------------------------------
### Plot filtering estimates
###-------------------------------------
df_sis <- tibble(
  t = 1:T,
  estimate = posterior_mean,
  lower = posterior_mean - 2 * posterior_sd,
  upper = posterior_mean + 2 * posterior_sd,
  theta_true = theta_true
)

plot_filtering_estimates <- function(df) {
  ggplot(df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Credible Interval"), alpha = 0.5) +
    geom_line(aes(y = theta_true, color = "True State"), linewidth = 0.5) +
    geom_line(aes(y = estimate, color = "Posterior Mean"), linewidth = 0.5) +
    scale_color_manual(values = c("True State" = "black", "Posterior Mean" = "blue")) +
    scale_fill_manual(values = c("95% Credible Interval" = "lightblue")) +
    ylab(TeX("$\\theta_t$")) +
    xlab("Time") +
    labs(color = "Line", fill = "Interval") +
    theme_bw() +
    theme(legend.position = "bottom")
}

# Display plot
plot_filtering_estimates(df_sis)

########################## Particle filter: Dynamic linear model ########################## 
# Load required packages
library(dplyr)
library(ggplot2)
library(latex2exp)

set.seed(10101)

# Parameters
n_particles <- 50000
sigma_w <- 1
sigma_mu <- 1
phi <- 0.5
T <- 200

###-----------------------------------
### Simulate true states and observations
###-----------------------------------
theta_true <- numeric(T)
y_obs <- numeric(T)

theta_true[1] <- rnorm(1, 0, sigma_w)
for (t in 2:T) {
  theta_true[t] <- rnorm(1, mean = theta_true[t - 1], sd = sigma_w)
}
y_obs <- rnorm(T, mean = phi * theta_true, sd = sigma_mu)

###-----------------------------------
### Initialize matrices for filtering
###-----------------------------------
particles <- matrix(0, n_particles, T)
particles_resampled <- matrix(0, n_particles, T)

weights <- matrix(0, n_particles, T)
weights_norm <- matrix(0, n_particles, T)
weights_fixed <- matrix(1 / n_particles, n_particles, T)

log_weights <- matrix(0, n_particles, T)

# Initialize t = 1
particles[, 1] <- rnorm(n_particles, 0, sigma_w)
weights[, 1] <- dnorm(y_obs[1], mean = phi * particles[, 1], sd = sigma_mu)
weights_norm[, 1] <- weights[, 1] / sum(weights[, 1])

# Resampling
resample_index <- sample(1:n_particles, size = n_particles, replace = TRUE, prob = weights_norm[, 1])
particles[, 1] <- particles[resample_index, 1]
particles_resampled[, 1] <- particles[, 1]

###-----------------------------------
### Particle filtering loop
###-----------------------------------
pb <- txtProgressBar(min = 1, max = T, style = 3)

for (t in 2:T) {
  # Propagate particles
  particles[, t] <- rnorm(n_particles, mean = particles[, t - 1], sd = sigma_w)
  
  # Log-weight and normalization
  log_weights[, t] <- dnorm(y_obs[t], mean = phi * particles[, t], sd = sigma_mu, log = TRUE)
  weights[, t] <- exp(log_weights[, t])
  weights_norm[, t] <- weights[, t] / sum(weights[, t])
  
  # Resampling
  resample_index <- sample(1:n_particles, size = n_particles, replace = TRUE, prob = weights_norm[, t])
  
  if (t < T) {
    particles[, 1:t] <- particles[resample_index, 1:t]
  } else {
    particles_resampled[, 1:t] <- particles[resample_index, 1:t]
  }
  
  setTxtProgressBar(pb, t)
}
close(pb)

###-----------------------------------
### Posterior summaries
###-----------------------------------
posterior_mean <- colSums(particles * weights_norm)
posterior_sd <- sqrt(colSums((particles^2) * weights_norm) - posterior_mean^2)

posterior_mean_resampled <- colSums(particles_resampled * weights_fixed)
posterior_sd_resampled <- sqrt(colSums((particles_resampled^2) * weights_fixed) - posterior_mean_resampled^2)

marginal_likelihood <- colMeans(weights)
plot(marginal_likelihood, type = "l", main = "Marginal Likelihood", xlab = "Time", ylab = "Likelihood")

###-----------------------------------
### Plot filtering results
###-----------------------------------
df_filter <- tibble(
  t = 1:T,
  mean = posterior_mean,
  lower = posterior_mean - 2 * posterior_sd,
  upper = posterior_mean + 2 * posterior_sd,
  mean_resampled = posterior_mean_resampled,
  lower_resampled = posterior_mean_resampled - 2 * posterior_sd_resampled,
  upper_resampled = posterior_mean_resampled + 2 * posterior_sd_resampled,
  theta_true = theta_true
)

plot_filtering_estimates <- function(df) {
  ggplot(df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "SIS 95% CI"), alpha = 0.3) +
    geom_line(aes(y = mean, color = "SIS Mean"), linewidth = 0.5) +
    geom_line(aes(y = mean_resampled, color = "Resampled Mean"), linewidth = 0.5) +
    geom_line(aes(y = theta_true, color = "True State"), linewidth = 0.5) +
    scale_color_manual(values = c("SIS Mean" = "blue", "Resampled Mean" = "purple", "True State" = "black")) +
    scale_fill_manual(values = c("SIS 95% CI" = "lightblue")) +
    labs(
      y = TeX("$\\theta_t$"), x = "Time",
      color = "Line", fill = "Band"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

plot_filtering_estimates(df_filter)

########################## Mining change point ########################## 

# Load libraries
library(coda)
library(latex2exp)

set.seed(10101)

# Load and inspect data
dataset <- read.csv("https://raw.githubusercontent.com/besmarter/BSTApp/refs/heads/master/DataApp/MiningDataCarlin.csv")
y <- dataset$Count
T_len <- length(y)

# Prior parameters
a10 <- 0.5; b10 <- 1
a20 <- 0.5; b20 <- 1

# MCMC settings
H <- 60
MCMC <- 20000
burnin <- 1000
S <- MCMC + burnin
keep <- (burnin + 1):S

# Preallocate
theta1 <- numeric(S)
theta2 <- numeric(S)
kk <- numeric(S)

pb <- txtProgressBar(min = 1, max = S, style = 3)

for (s in 1:S) {
  a1 <- a10 + sum(y[1:H])
  b1 <- b10 + H
  theta1[s] <- rgamma(1, a1, b1)
  
  a2 <- a20 + sum(y[(H + 1):T_len])
  b2 <- b20 + T_len - H
  theta2[s] <- rgamma(1, a2, b2)
  
  log_p <- numeric(T_len)
  for (l in 1:T_len) {
    log_p[l] <- l * (theta2[s] - theta1[s]) + sum(y[1:l]) * log(theta1[s] / theta2[s])
  }
  p <- exp(log_p - max(log_p))  # Stabilize numerics
  prob <- p / sum(p)
  H <- sample(1:T_len, 1, prob = prob)
  kk[s] <- H
  
  setTxtProgressBar(pb, s)
}
close(pb)

# Post-burnin chains
theta1_post <- mcmc(theta1[keep])
theta2_post <- mcmc(theta2[keep])
H_post <- mcmc(kk[keep])

# Diagnostics
summary(theta2_post)
plot(theta2_post, density = FALSE, main = "Trace plot", ylab = TeX("$\\theta_{2}$"))
autocorr.plot(theta2_post)
raftery.diag(theta2_post)
geweke.diag(theta2_post)
heidel.diag(theta2_post)
effectiveSize(theta2_post)

# -------------------------------
# Marginal-conditional simulator
# -------------------------------
theta1_prior <- rgamma(MCMC, a10, b10)
theta2_prior <- rgamma(MCMC, a20, b20)
k_prior <- sample(1:T_len, MCMC, replace = TRUE)

simulate_y <- function(par) {
  y1 <- rpois(par[3], par[1])
  y2 <- if (par[3] < T_len) rpois(T_len - par[3], par[2]) else numeric(0)
  c(y1, y2)
}

pars_mc <- cbind(theta1_prior, theta2_prior, k_prior)
Yt <- apply(pars_mc, 1, simulate_y)
mcmc_mc <- mcmc(pars_mc)
summary_mc <- summary(mcmc_mc)

# -------------------------------
# Successive-conditional simulator
# -------------------------------
simulate_succ <- function(a10, b10, a20, b20, par) {
  y <- simulate_y(par)
  H <- par[3]
  
  a1 <- a10 + sum(y[1:H])
  b1 <- b10 + H
  theta1_new <- rgamma(1, a1, b1)
  
  if (H == T_len) {
    a2 <- a20
  } else {
    a2 <- a20 + sum(y[(H + 1):T_len])
  }
  b2 <- b20 + T_len - H
  theta2_new <- rgamma(1, a2, b2)
  
  log_p <- sapply(1:T_len, function(l) {
    l * (theta2_new - theta1_new) + sum(y[1:l]) * log(theta1_new / theta2_new)
  })
  p <- exp(log_p - max(log_p))
  prob <- p / sum(p)
  H_new <- sample(1:T_len, 1, prob = prob)
  
  list(y = y, pars = c(theta1_new, theta2_new, H_new))
}

# Run simulation
a10 <- 0.5; b10 <- 1
a20 <- 0.5; b20 <- 1
pars_sc <- matrix(NA, nrow = MCMC, ncol = 3)
pars_sc[1, ] <- c(rgamma(1, a10, b10), rgamma(1, a20, b20), sample(1:T_len, 1))

for (s in 2:MCMC) {
  res <- simulate_succ(a10, b10, a20, b20, pars_sc[s - 1, ])
  pars_sc[s, ] <- res$pars
}

mcmc_sc <- mcmc(pars_sc)
summary_sc <- summary(mcmc_sc)

# -------------------------------
# Geweke test
# -------------------------------
geweke_test <- function(j) {
  num <- summary_mc$statistics[j, 1] - summary_sc$statistics[j, 1]
  denom <- sqrt(summary_mc$statistics[j, 4]^2 + summary_sc$statistics[j, 4]^2)
  z <- num / denom
  reject <- abs(z) > qnorm(0.975)
  list(Test = z, Reject = reject)
}

geweke_test(1); geweke_test(2); geweke_test(3)

# -------------------------------
# Repeat with wrong prior to test sensitivity
# -------------------------------
a10 <- 1; b10 <- 0.5
a20 <- 1; b20 <- 0.5
pars_sc2 <- matrix(NA, nrow = MCMC, ncol = 3)
pars_sc2[1, ] <- c(rgamma(1, a10, b10), rgamma(1, a20, b20), sample(1:T_len, 1))

for (s in 2:MCMC) {
  res <- simulate_succ(a10, b10, a20, b20, pars_sc2[s - 1, ])
  pars_sc2[s, ] <- res$pars
}

mcmc_sc2 <- mcmc(pars_sc2)
summary_sc2 <- summary(mcmc_sc2)

geweke_test <- function(j) {
  num <- summary_mc$statistics[j, 1] - summary_sc2$statistics[j, 1]
  denom <- sqrt(summary_mc$statistics[j, 4]^2 + summary_sc2$statistics[j, 4]^2)
  z <- num / denom
  reject <- abs(z) > qnorm(0.975)
  list(Test = z, Reject = reject)
}

geweke_test(1); geweke_test(2); geweke_test(3)

