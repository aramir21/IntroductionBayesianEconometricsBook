#### 401k: Treatment effects ####
rm(list = ls())
set.seed(10101)

library(MCMCpack)
library(coda)
library(ggplot2)

mydata <- read.csv("https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/401k.csv", sep = ",", header = TRUE, quote = "")
attach(mydata )
y <- net_tfa 
# net_tfa: net financial assets
# Regressors quantity including intercept
X <- cbind(e401, age, inc, fsize, educ, marr, twoearn, db, pira, hown, 1)
# e401: 401k eligibility 
# age, income, family size, years of education, marital status indicator
# two-earner status indicator, defined benefit pension status indicator
# IRA participation indicator, and home ownership indicator.

# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
posterior  <- MCMCpack::MCMCregress(y~X-1)
summary(coda::mcmc(posterior))

# Extract posterior samples for e401 (first coefficient)
beta_e401 <- posterior[, 1]

# Convert to data frame for ggplot
df_posterior <- data.frame(beta_e401)

# Plot with ggplot
ggplot(df_posterior, aes(x = beta_e401)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = mean(beta_e401), color = "red", linetype = "dashed") +
  labs(
    title = "Posterior Distribution of Treatment Effect (e401)",
    x = expression(beta[e401]),
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

##### Collider bias ######
rm(list = ls())
set.seed(10101)

library(MCMCpack)
library(dplyr)

# Parameters
N <- 5000        # Sample size
R <- 100         # Number of replications
true_effect <- 2 # True causal effect
burn <- 1000
mcmc_draws <- 3000

# Storage for posterior means
results <- data.frame(
  rep = 1:R,
  correct = numeric(R),
  biased = numeric(R)
)

# Simulation loop
for (r in 1:R) {
  # Simulate data based on DAG
  X <- rnorm(N, 0, 1)
  D <- 1 + 0.5 * X + rnorm(N, 0, 0.1)
  C <- 1 + 1 * D + 1 * X + rnorm(N, 0, 0.1) > 2
  Y <- 2 + true_effect * D + 0.1 * X + rnorm(N, 0, 0.1) # Strong X->Y path
  
  # Fit Bayesian regressions
  model_correct <- MCMCregress(Y ~ D + X, burnin = burn, mcmc = mcmc_draws, verbose = FALSE)
  model_biased  <- MCMCregress(Y ~ D + X + C, burnin = burn, mcmc = mcmc_draws, verbose = FALSE)
  
  # Extract posterior means for D
  results$correct[r] <- mean(as.matrix(model_correct)[, "D"])
  results$biased[r]  <- mean(as.matrix(model_biased)[, "D"])
}

# Compute bias for each model
results <- results %>%
  mutate(
    bias_correct = correct - true_effect,
    bias_biased = biased - true_effect
  )

# Average bias
avg_bias <- results %>%
  summarise(
    mean_correct = mean(bias_correct),
    mean_biased = mean(bias_biased),
    sd_correct = sd(bias_correct),
    sd_biased = sd(bias_biased)
  )

print(avg_bias)

# Optional: Visualize distribution of estimates across replications
library(ggplot2)
df_long <- results %>%
  select(rep, correct, biased) %>%
  tidyr::pivot_longer(cols = c(correct, biased), names_to = "model", values_to = "estimate")

ggplot(df_long, aes(x = estimate, fill = model)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = true_effect, color = "black", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Distribution of Posterior Means Across 100 Simulations",
    x = expression(paste("Posterior Mean of ", beta[D])),
    y = "Density"
  ) +
  theme_minimal(base_size = 14)
