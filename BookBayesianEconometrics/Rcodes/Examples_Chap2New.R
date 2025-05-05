# Examples in chapter 2


########################## The likelihood principle ########################## 
# Binomial test: one-sided significance level for observing 9 or more successes in 12 trials

# Parameters
successes <- 9         # Number of observed successes
n_trials <- 12         # Total number of trials
p_null <- 0.5          # Null hypothesis success probability

# Calculate one-sided significance level
significance_level <- sum(dbinom(successes:n_trials, size = n_trials, prob = p_null))

# Output result with context
message(sprintf("One-sided significance level (P(X ≥ %d | H0: p = %.1f)): %.4f", 
                successes, p_null, significance_level))

# Negative binomial test: Probability of observing at least 3 tails before 9 heads (failures)

# Parameters
target_successes <- 3      # Number of target successes (e.g., tails)
num_failures <- 9          # Number of failures (e.g., heads)
p_success <- 0.5           # Probability of success (e.g., tails)

# Compute the one-sided significance level: P(X ≥ 3 successes before 9 failures)
significance_level <- 1 - pnbinom(q = num_failures - 1,
                                  size = target_successes,
                                  prob = p_success)

# Print result
message(sprintf("P(at least %d tails before %d heads): %.4f",
                target_successes, num_failures, significance_level))

########################## A simple working example ##########################

# Load required package
library(metRology)

# Input values
n <- 50                # Sample size
sample_mean <- 102     # Sample mean
sample_var <- 10       # Sample variance
alpha <- n - 1         # Degrees of freedom
se <- sqrt(sample_var / n)  # Standard error

# Confidence interval for the mean (95%)
lower_bound <- sample_mean - abs(qt(0.025, df = alpha)) * se
upper_bound <- sample_mean + abs(qt(0.025, df = alpha)) * se

cat(sprintf("95%% CI: [%.3f, %.3f]\n", lower_bound, upper_bound))

# Probability that population mean is greater than y_cut
y_cut <- 103
prob_mu_gt_ycut <- 1 - pt.scaled(y_cut, df = alpha, mean = sample_mean, sd = se)

cat(sprintf("P(mu > %.1f) = %.4f\n", y_cut, prob_mu_gt_ycut))

