# Examples in chapter 1


########################## The Base rate fallacy ########################## 

# Define known probabilities
prob_disease <- 0.002        # P(Disease)
sensitivity <- 0.98          # P(Positive test | Disease) - True Positive Rate
false_positive_rate <- 0.02  # P(Positive test | No disease) - False Positive Rate

# Compute posterior using Bayes' theorem
posterior_prob <- (prob_disease * sensitivity) /
  (prob_disease * sensitivity + (1 - prob_disease) * false_positive_rate)

# Print the result
message(sprintf("Probability of disease given a positive test is %.2f", posterior_prob))

##########################  Jesus Crist ##########################

# Define known probabilities
prior_resurrection <- 2 / (117 * 1e9)  # P(Resurrection)
true_positive_rate <- 0.9999999         # P(Witness | Resurrection)

# Assume a 50% false positive rate as implied in earlier context
false_positive_rate <- 0.5              # P(Witness | No Resurrection)

# Compute posterior probability using Bayes' rule
posterior_resurrection <- (prior_resurrection * true_positive_rate) /
  (prior_resurrection * true_positive_rate +
     (1 - prior_resurrection) * false_positive_rate)

# Print result
message(sprintf("Probability of resurrection given a witness is %.2e", posterior_resurrection))

##########################  The Monty Hall problem ########################## 
# Set simulation seed for reproducibility
set.seed(010101)

# Number of simulations
num_simulations <- 100000

# Monty Hall game function
simulate_game <- function(switch_door = FALSE) {
  doors <- 1:3
  car_location <- sample(doors, 1)
  first_guess <- sample(doors, 1)
  
  # Host reveals a goat
  if (car_location != first_guess) {
    host_opens <- doors[!doors %in% c(car_location, first_guess)]
  } else {
    host_opens <- sample(doors[doors != first_guess], 1)
  }
  
  # Determine second guess if player switches
  second_guess <- doors[!doors %in% c(first_guess, host_opens)]
  
  win_if_no_switch <- (first_guess == car_location)
  win_if_switch <- (second_guess == car_location)
  
  if (switch_door) {
    return(win_if_switch)
  } else {
    return(win_if_no_switch)
  }
}

# Simulate without switching
prob_no_switch <- mean(replicate(num_simulations, simulate_game(switch_door = FALSE)))
message(sprintf("Winning probability without switching: %.3f", prob_no_switch))

# Simulate with switching
prob_with_switch <- mean(replicate(num_simulations, simulate_game(switch_door = TRUE)))
message(sprintf("Winning probability with switching: %.3f", prob_with_switch))

# Poisson-Gamma Bayesian Inference with Vague and Empirical Bayes Priors

# Load required packages
library(ggplot2)
library(ggpubr)
library(latex2exp)

# Set seed and data
y <- c(0, 3, 2, 1, 0)
n_obs <- length(y)
set.seed(10101)

# Summary statistics
message(sprintf("Sample mean: %.2f", mean(y)))
message(sprintf("Sample variance: %.2f", var(y)))

# Predictive probability function
predictive_prob <- function(y, a0, b0) {
  n <- length(y)
  an <- a0 + sum(y)
  bn <- b0 / (b0 * n + 1)
  p <- bn / (bn + 1)
  prob <- 1 - pnbinom(0, size = an, prob = 1 - p)
  return(prob)
}

# Vague prior parameters
a0_vague <- 0.001
b0_vague <- 1 / a0_vague

prior_mean_vague <- a0_vague * b0_vague
prior_var_vague  <- a0_vague * b0_vague^2

message(sprintf("Vague prior mean: %.2f, variance: %.2f", prior_mean_vague, prior_var_vague))

# Predictive under vague prior
predictive_vague <- predictive_prob(y, a0 = a0_vague, b0 = b0_vague)
message(sprintf("P(y > 0) under vague prior: %.2f", predictive_vague))

# Empirical Bayes: log marginal likelihood
log_marginal_likelihood <- function(theta, y) {
  a0 <- theta[1]
  b0 <- theta[2]
  n <- length(y)
  if (a0 <= 0 || b0 <= 0) return(Inf)
  an <- a0 + sum(y)
  loglik <- lgamma(an) + sum(y) * log(b0 / (n * b0 + 1)) - a0 * log(n * b0 + 1) - lgamma(a0)
  return(-loglik)
}

# Optimize
theta_start <- c(0.01, 1 / 0.1)
opt_result <- optim(
  par = theta_start,
  fn = log_marginal_likelihood,
  method = "BFGS",
  control = list(maxit = 1000),
  hessian = TRUE,
  y = y
)

# Extract parameters
a0_eb <- opt_result$par[1]
b0_eb <- opt_result$par[2]

prior_mean_eb <- a0_eb * b0_eb
prior_var_eb  <- a0_eb * b0_eb^2

message(sprintf("Empirical Bayes prior: shape = %.4f, scale = %.4f", a0_eb, b0_eb))
message(sprintf("Empirical Bayes mean: %.4f, variance: %.4f", prior_mean_eb, prior_var_eb))

# Predictive under EB prior
predictive_eb <- predictive_prob(y, a0 = a0_eb, b0 = b0_eb)
message(sprintf("P(y > 0) under EB prior: %.4f", predictive_eb))

# Plot densities
lambda_vals <- seq(0.01, 10, 0.01)

prior_vague <- dgamma(lambda_vals, shape = a0_vague, scale = b0_vague)
prior_eb <- dgamma(lambda_vals, shape = a0_eb, scale = b0_eb)
posterior_vague <- dgamma(lambda_vals, shape = a0_vague + sum(y), scale = b0_vague / (b0_vague * n_obs + 1))
posterior_eb <- dgamma(lambda_vals, shape = a0_eb + sum(y), scale = b0_eb / (b0_eb * n_obs + 1))

# Likelihood function
likelihood_vals <- sapply(lambda_vals, function(lam) prod(dpois(y, lam)))
likelihood_scaled <- likelihood_vals * max(posterior_eb) / max(likelihood_vals)

# Prepare data
densities_df <- data.frame(
  lambda = lambda_vals,
  VaguePrior = prior_vague,
  EBPrior = prior_eb,
  PosteriorV = posterior_vague,
  PosteriorEB = posterior_eb,
  Likelihood = likelihood_scaled
)

# Plots
plot_density <- function(data, yval, title) {
  ggplot(data, aes(x = lambda, y = !!as.name(yval))) +
    geom_line() +
    xlab(TeX("$lambda$")) +
    ylab("Density") +
    ggtitle(title)
}

figs <- list(
  plot_density(densities_df, "VaguePrior", "Prior: Vague Gamma"),
  plot_density(densities_df, "EBPrior", "Prior: Empirical Bayes Gamma"),
  plot_density(densities_df, "PosteriorV", "Posterior: Vague Gamma"),
  plot_density(densities_df, "PosteriorEB", "Posterior: Empirical Bayes Gamma")
)

combined_plot <- ggarrange(plotlist = figs, ncol = 2, nrow = 2)
annotate_figure(combined_plot, top = text_grob("Vague vs Empirical Bayes: Poisson-Gamma", face = "bold", size = 14))

# Prior, Likelihood, and Posterior comparison (Empirical Bayes)
dataNew <- data.frame(cbind(rep(lambda_vals, 3), c(prior_eb, posterior_eb, likelihood_scaled),
                            rep(1:3, each = length(lambda_vals))))
colnames(dataNew) <- c("Lambda", "Density", "Factor")
dataNew$Factor <- factor(dataNew$Factor, levels = c("1", "3", "2"),
                         labels = c("Prior", "Likelihood", "Posterior"))

ggplot(data = dataNew, aes(x = Lambda, y = Density, group = Factor, color = Factor)) +
  geom_line() +
  xlab(TeX("$lambda$")) +
  ylab("Density") +
  ggtitle("Prior, Likelihood, and Posterior: Empirical Bayes Poisson-Gamma Model") +
  guides(color = guide_legend(title = "Information")) +
  scale_color_manual(values = c("red", "yellow", "blue"))

# Predictive densities
predictive_density <- function(y, y0, a0, b0) {
  n <- length(y)
  an <- a0 + sum(y)
  bn <- b0 / (b0 * n + 1)
  p <- bn / (bn + 1)
  dnbinom(y0, size = an, prob = 1 - p)
}

y0_vals <- 0:10
pred_vague <- predictive_density(y, y0_vals, a0_vague, b0_vague)
pred_eb <- predictive_density(y, y0_vals, a0_eb, b0_eb)

pred_df <- data.frame(y0 = y0_vals, Vague = pred_vague, EmpBayes = pred_eb)

ggplot(pred_df) +
  geom_point(aes(x = y0, y = Vague, color = "Vague")) +
  geom_point(aes(x = y0, y = EmpBayes, color = "Empirical Bayes")) +
  xlab(TeX("$y_0$")) +
  ylab("Density") +
  ggtitle("Predictive Densities: Vague vs Empirical Bayes Priors") +
  scale_color_manual(name = "Prior", values = c("Vague" = "red", "Empirical Bayes" = "blue")) +
  theme_minimal()

# Posterior odds and model probabilities
log_mg_vague <- -log_marginal_likelihood(c(a0_vague, b0_vague), y)
log_mg_eb <- -log_marginal_likelihood(c(a0_eb, b0_eb), y)

odds_eb_vs_vague <- exp(log_mg_eb - log_mg_vague)
prob_eb <- odds_eb_vs_vague / (1 + odds_eb_vs_vague)
prob_vague <- 1 - prob_eb

message(sprintf("Posterior model probabilities - Empirical Bayes: %.4f, Vague: %.4f", prob_eb, prob_vague))

# BMA: mean and variance
post_mean_eb <- (a0_eb + sum(y)) * (b0_eb / (b0_eb * n_obs + 1))
post_var_eb <- (a0_eb + sum(y)) * (b0_eb / (b0_eb * n_obs + 1))^2
post_mean_vague <- (a0_vague + sum(y)) * (b0_vague / (b0_vague * n_obs + 1))
post_var_vague <- (a0_vague + sum(y)) * (b0_vague / (b0_vague * n_obs + 1))^2

bma_mean <- prob_eb * post_mean_eb + prob_vague * post_mean_vague
bma_var <- prob_eb * post_var_eb + prob_vague * post_var_vague +
  prob_eb * (post_mean_eb - bma_mean)^2 + prob_vague * (post_mean_vague - bma_mean)^2

message(sprintf("BMA posterior mean: %.4f, variance: %.4f", bma_mean, bma_var))

# BMA Predictive
pred_bma <- prob_eb * pred_eb + prob_vague * pred_vague
pred_df_bma <- data.frame(y0 = y0_vals, BMAPredictive = pred_bma)

ggplot(data = pred_df_bma) + 
  geom_point(aes(y0, BMAPredictive, color = "red")) +  
  xlab(TeX("$y_0$")) + 
  ylab("Density") + 
  ggtitle("Predictive density: BMA") + 
  guides(color = guide_legend(title = "BMA")) + 
  scale_color_manual(labels = c("Probability"), values = c("red")) + 
  scale_x_continuous(breaks = seq(0, 10, by = 1))

# Sequential Bayesian updating
sequential_update <- function(y, lambda, a0, b0) {
  n <- length(y)
  results <- matrix(NA, nrow = length(lambda), ncol = n)
  a_new <- a0
  b_new <- b0
  
  for (i in 1:n) {
    a_new <- a_new + y[i]
    b_new <- b_new / (b_new * 1 + 1)
    results[, i] <- dgamma(lambda, shape = a_new, scale = b_new)
  }
  return(results)
}

lambda_grid <- lambda_vals
updates <- sequential_update(y, lambda_grid, a0_vague, b0_vague)

updates_df <- data.frame(
  lambda = rep(lambda_grid, times = n_obs),
  density = as.vector(updates),
  iteration = factor(rep(1:n_obs, each = length(lambda_grid)))
)

ggplot(updates_df, aes(x = lambda, y = density, color = iteration)) +
  geom_line() +
  xlab(TeX("$lambda$")) +
  ylab("Density") +
  ggtitle("Sequential Bayesian Updating: Vague Prior") +
  theme_minimal() +
  guides(color = guide_legend(title = "Update step"))

########################## Decision theory ########################## 

# Load required package
library(HDInterval)

# Posterior parameters
posterior_shape <- sum(y) + a0_eb
posterior_scale <- b0_eb / (n_obs * b0_eb + 1)

# Posterior draws
n_draws <- 1e6
posterior_draws <- rgamma(n_draws, shape = posterior_shape, scale = posterior_scale)

# Optimal point estimates
mean_quad_loss <- posterior_shape * posterior_scale
median_abs_loss <- qgamma(0.5, shape = posterior_shape, scale = posterior_scale)

# Generalized absolute loss with asymmetric costs
k0 <- 2
k1 <- 1
quantile_weight <- k0 / (k0 + k1)
generalized_abs_loss <- quantile(posterior_draws, probs = quantile_weight)

# Display point estimates
message(sprintf("Quadratic loss (mean): %.4f", mean_quad_loss))
message(sprintf("Absolute loss (median): %.4f", median_abs_loss))
message(sprintf("Generalized absolute loss (quantile %.2f): %.4f", quantile_weight, generalized_abs_loss))

# Hypothesis test: H0: λ ∈ [0, 1), H1: λ ≥ 1
k0 <- 1
k1 <- 1
prob_h0 <- pgamma(1, shape = posterior_shape, scale = posterior_scale)
prob_h1 <- 1 - prob_h0
decision_threshold <- k1 / (k0 + k1)

message(sprintf("P(H0: λ < 1): %.4f", prob_h0))
message(sprintf("P(H1: λ ≥ 1): %.4f", prob_h1))
message(sprintf("Reject H0? %s", ifelse(prob_h1 > decision_threshold, "Yes", "No")))

# Credible intervals
ci_lower <- qgamma(0.025, shape = posterior_shape, scale = posterior_scale)
ci_upper <- qgamma(0.975, shape = posterior_shape, scale = posterior_scale)
hdi_interval <- HDInterval::hdi(posterior_draws, credMass = 0.95)

message(sprintf("95%% Credible Interval (Quantiles): [%.4f, %.4f]", ci_lower, ci_upper))
message(sprintf("95%% HDI: [%.4f, %.4f]", hdi_interval[1], hdi_interval[2]))

# Optimal prediction under quadratic loss
nbinom_p <- posterior_scale / (posterior_scale + 1)
opt_predictive <- (nbinom_p / (1 - nbinom_p)) * posterior_shape
message(sprintf("Optimal predictive (quadratic loss): %.4f", opt_predictive))

