#### Basic Example ####
set.seed(10101)

# Parameters
mu1 <- 15; mu0 <- 10
sigma1 <- sqrt(4); sigma0 <- sqrt(2)
N1 <- 100; N0 <- 100

# Simulate data
y1 <- rnorm(N1, mu1, sigma1)
y0 <- rnorm(N0, mu0, sigma0)

# Prior hyperparameters
alpha0 <- 0.01
delta0 <- 0.01
beta0 <- 0  # Not used now

simulate_t <- function(y) {
  N <- length(y)
  ybar <- mean(y)
  sse <- sum((y - ybar)^2)
  
  alpha_n <- alpha0 + N
  delta_n <- sse + delta0
  scale2 <- ((N + 1) * delta_n) / (N * alpha_n)
  
  df <- alpha_n
  loc <- ybar
  scale <- sqrt(scale2)
  
  rt(1, df = df) * scale + loc
}

# Posterior predictive draws
ppd1 <- replicate(1000, simulate_t(y1))
ppd0 <- replicate(1000, simulate_t(y0))

# Summaries
mean(ppd1); mean(ppd0)

# Plot
hist(ppd1, col = rgb(1, 0, 0, 0.4), freq = FALSE, main = "Posterior Predictive",
     xlab = "Y", xlim = c(5, 25))
hist(ppd0, col = rgb(0, 0, 1, 0.4), freq = FALSE, add = TRUE)
legend("topright", legend = c("Treatment", "Control"),
       fill = c(rgb(1, 0, 0, 0.4), rgb(0, 0, 1, 0.4)))

# Posterior distribution of ATE
posterior_mu <- function(y) {
  N <- length(y)
  ybar <- mean(y)
  sse <- sum((y - ybar)^2)
  
  alpha_n <- alpha0 + N
  delta_n <- sse + delta0
  
  # Posterior variance for mean
  var_mu <- delta_n / (alpha_n * N)
  df <- alpha_n
  
  # Draw from t for parameter μ
  loc <- ybar
  scale <- sqrt(var_mu)
  rt(1, df = df) * scale + loc
}

# Posterior draws for parameters
n_draws <- 10000
mu1_draws <- replicate(n_draws, posterior_mu(y1))
mu0_draws <- replicate(n_draws, posterior_mu(y0))

# Parameter uncertainty: difference of means
ate_draws <- mu1_draws - mu0_draws
summary(coda::mcmc(ate_draws))
# Summaries
ate_mean <- mean(ate_draws)
ci <- quantile(ate_draws, c(0.025, 0.975))

cat("Posterior mean of ATE:", ate_mean, "\n")
cat("95% Credible Interval:", ci, "\n")

# Plot posterior distribution of ATE
hist(ate_draws, breaks = 50, freq = FALSE,
     main = "Posterior Distribution of ATE",
     xlab = "ATE (Y(1) - Y(0))", col = "lightblue", border = "white")
abline(v = ate_mean, col = "red", lwd = 2)
abline(v = ci, col = "darkgreen", lty = 2, lwd = 2)

legend("topright", legend = c("Posterior Mean", "95% Credible Interval"),
       col = c("red", "darkgreen"), lwd = 2, lty = c(1, 2),
       bty = "n", cex = 0.8)  # Smaller legend using cex

y <- c(y0, y1)
d <- c(rep(0, N0), rep(1, N1))

data_reg <- data.frame(y = y, d = d)

# OLS regression
model <- lm(y ~ d, data = data_reg)
summary(model)

#### 401k: Treatment effects ####
rm(list = ls())
set.seed(10101)

library(MCMCpack)
library(coda)
library(ggplot2)

mydata <- read.csv("https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/401k.csv", sep = ",", header = TRUE, quote = "")
attach(mydata )
y <- net_tfa/1000 
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

#### Collider bias ####
rm(list = ls())
set.seed(10101)

library(MCMCpack)
library(dplyr)
library(ggplot2)

# Parameters
n <- 1000
true_effect <- 2
replications <- 100

# Store results
results <- data.frame(
  rep = 1:replications,
  correct = numeric(replications),
  biased = numeric(replications)
)

for (r in 1:replications) {
  
  # Simulate data under DAG: S -> B -> Y, S -> Y, U -> B, U -> Y
  U <- rnorm(n, 0, 1)
  S <- rbinom(n, prob = 0.7, size = 1)
  B <- S + U + rnorm(n)
  Y <- S + B + 1.5*U + rnorm(n)
  
  # Correct model: does NOT condition on collider B
  model_correct <- MCMCregress(Y ~ S, burnin = 1000, mcmc = 3000, verbose = FALSE)
  
  # Biased model: conditions on collider B
  model_biased <- MCMCregress(Y ~ S + B, burnin = 1000, mcmc = 3000, verbose = FALSE)
  
  # Posterior means for S
  results$correct[r] <- mean(as.matrix(model_correct)[, "S"])
  results$biased[r] <- mean(as.matrix(model_biased)[, "S"])
}

# Compute bias
results <- results %>%
  mutate(
    bias_correct = correct - true_effect,
    bias_biased = biased - true_effect
  )

# Average bias and SD
avg_bias <- results %>%
  summarise(
    mean_correct = mean(bias_correct),
    mean_biased = mean(bias_biased),
    sd_correct = sd(bias_correct),
    sd_biased = sd(bias_biased)
  )

print(avg_bias)

# Visualization: distribution of posterior means across 100 simulations
df_long <- results %>%
  select(rep, correct, biased) %>%
  tidyr::pivot_longer(cols = c(correct, biased), names_to = "model", values_to = "estimate")

ggplot(df_long, aes(x = estimate, fill = model)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = true_effect, color = "black", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Posterior Means Across 100 Simulations",
    x = expression(paste("Posterior Mean of ", beta[S])),
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

######## Instrumental variables: 401k participation #############
# VERY IMPORTANT: The rivGibbs function from bayesm does not control for
#                 exogenous regressors in the first stage equation.
#                 This implies that the Wald ratio is not well calculated.
#                 This is way I perform from scratch the Gibbs sampler!!!
rm(list = ls())
set.seed(10101)
# Load packages
library(coda)
library(ggplot2)

# Load data
mydata <- read.csv("https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/401k.csv",
                   sep = ",", header = TRUE, quote = "")
# Attach variables
attach(mydata)
y <- net_tfa/1000                  # Outcome: net financial assets
x <- as.vector(p401)          # Endogenous regressor: participation
w <- as.matrix(cbind(1, age, inc, fsize, educ, marr, twoearn, db, pira, hown))  # Exogenous regressors with intercept
z <- as.matrix(e401)          # Instrument: eligibility (NO intercept here)
X <- cbind(x, w)
Z <- cbind(z, w)
# Dimensions
k <- ncol(X)
kz <- ncol(Z)  
# Priors
b0 <- rep(0, k)
B0i <- diag(1e-5, k)
g0 <- rep(0, kz)
G0i <- diag(1e-5, kz)
nu <- 3
Psi0 <- nu * 1000 * diag(2)
Psi0i <- solve(Psi0)
# MCMC parameters
mcmc <- 5000
burnin <- 1000
tot <- mcmc + burnin
thin <- 1
# Auxiliary elements
XtX <- t(X)%*%X 
ZtZ <- t(Z)%*%Z 
nun <- nu + length(y)
# Gibbs sampling
PostBeta <- function(Sigma, Gamma){
  w1 <- Sigma[1,1] - Sigma[1,2]^2/Sigma[2,2]
  Bn <- solve(w1^(-1)*XtX + B0i)
  yaux <- y - (Sigma[1,2]/Sigma[2,2])*(x - Z%*%Gamma)
  bn <- Bn%*%(B0i%*%b0 + w1^(-1)*t(X)%*%yaux)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
PostGamma <- function(Sigma, Beta){
  w2 <- Sigma[2,2] - Sigma[1,2]^2/Sigma[1,1]
  Gn <- solve(w2^(-1)*ZtZ + G0i)
  xaux <- x - (Sigma[1,2]/Sigma[1,1])*(y - X%*%Beta)
  gn <- Gn%*%(G0i%*%g0 + w2^(-1)*t(Z)%*%xaux)
  Gamma <- MASS::mvrnorm(1, gn, Gn)
  return(Gamma)
}
PostSigma <- function(Beta, Gamma){
  Uy <- y - X%*%Beta; Ux <- x - Z%*%Gamma
  U <- cbind(Uy, Ux)
  Psin <- solve(Psi0i + t(U)%*%U)
  Sigmai <- rWishart::rWishart(1, df = nun, Sigma = Psin)
  Sigma <- solve(Sigmai[,,1]) 
  return(Sigma)
}
PostBetas <- matrix(0, tot, k)
PostGammas <- matrix(0, tot, kz)
PostSigmas <- matrix(0, tot, 2*(2+1)/2)
Beta <- rep(0, k); Gamma <- rep(0, kz)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Sigma <- PostSigma(Beta = Beta, Gamma = Gamma)
  Beta <- PostBeta(Sigma = Sigma, Gamma = Gamma)
  Gamma <- PostGamma(Sigma = Sigma, Beta = Beta)
  PostBetas[s,] <- Beta
  PostGammas[s,] <- Gamma
  PostSigmas[s,] <- matrixcalc::vech(Sigma)
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
Bs <- PostBetas[keep,]
Gs <- PostGammas[keep,]
Sigmas <- PostSigmas[keep,]
summary(coda::mcmc(Bs))
summary(coda::mcmc(Gs))
summary(coda::mcmc(Sigmas))

# Extract posterior draws for the treatment effect (participation = p401)
beta_draws <- Bs[,1]

# Convert to data frame for ggplot
df_beta <- data.frame(effect = as.vector(beta_draws))

# Plot posterior distribution of treatment effect
ggplot(df_beta, aes(x = effect)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = mean(beta_draws), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Posterior Distribution of 401(k) Participation Effect",
    x = expression(beta["p401"]),
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

ITT <- Bs[,1]*Gs[,1]
summary(coda::mcmc(ITT))
# Convert to data frame for ggplot
df_ITT <- data.frame(ITT = as.vector(ITT))

# Plot posterior distribution of treatment effect
ggplot(df_ITT, aes(x = ITT)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = mean(ITT), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Posterior Distribution of 401(k) ITT",
    x = expression(beta["ITT"]),
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

# library(AER)
# fs <- lm(x ~ z - 1)
# summary(fs)
# iv_model <- ivreg(formula = y ~ x + w | w + z, data = mydata)
# summary(iv_model, diagnostics = TRUE)
# 
# rf <- lm(net_tfa ~ e401 + age + inc + fsize + educ + marr + twoearn + db + pira + hown, data = mydata)
# fs <- lm(p401 ~ e401 + age + inc + fsize + educ + marr + twoearn + db + pira + hown, data = mydata)
# 
# conditional_wald <- coef(rf)["e401"] / coef(fs)["e401"]
# conditional_wald
# # Controls
# controls <- c("age", "inc", "fsize", "educ", "marr", "twoearn", "db", "pira", "hown")
# 
# # Step 1: Partial out controls
# Y_res <- lm(net_tfa ~ ., data = mydata[, controls])$residuals
# D_res <- lm(p401 ~ ., data = mydata[, controls])$residuals
# Z_res <- lm(e401 ~ ., data = mydata[, controls])$residuals
# 
# # Step 2: Conditional Wald ratio
# wald_partial <- cov(Y_res, Z_res) / cov(D_res, Z_res)
# wald_partial

