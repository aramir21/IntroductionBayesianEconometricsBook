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
     xlab = "Y(d)", xlim = c(5, 25))
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
     xlab = "ATE", col = "lightblue", border = "white")
abline(v = ate_mean, col = "red", lwd = 2)
abline(v = ci, col = "darkgreen", lty = 2, lwd = 2)

legend("topright", legend = c("Posterior Mean", "95% Credible Interval"),
       col = c("red", "darkgreen"), lwd = 2, lty = c(1, 2),
       bty = "n", cex = 0.8)  # Smaller legend using cex

y <- c(y0, y1)
d <- c(rep(0, N0), rep(1, N1))

data_reg <- data.frame(y = y, d = d)

# # OLS regression
# model <- lm(y ~ d, data = data_reg)
# summary(model)

#### Vitamin A supplements #####
rm(list = ls())
set.seed(10101)
library(dplyr)
# Simulate data
Nc111 <- 9663
c111 <- cbind(rep(1, Nc111), rep(1, Nc111), rep(1, Nc111)) 
Nc110 <- 12
c110 <- cbind(rep(1, Nc110), rep(1, Nc110), rep(0, Nc110)) 
Nn101 <- 2385
n101 <- cbind(rep(1, Nn101), rep(0, Nn101), rep(1, Nn101)) 
Nn100 <- 34
n100 <- cbind(rep(1, Nn100), rep(0, Nn100), rep(0, Nn100)) 
Ncn001 <- 11514
cn001 <- cbind(rep(0, Ncn001), rep(0, Ncn001), rep(1, Ncn001)) 
Ncn000 <- 74
cn000 <- cbind(rep(0, Ncn000), rep(0, Ncn000), rep(0, Ncn000)) 

mydata <- rbind(c111, c110, n101, n100, cn001, cn000)
mydata <- data.frame(Z = mydata[,1], D = mydata[,2], Y = mydata[,3])
N <- dim(mydata)[1]
attach(mydata)

# Sampling function C (type)
SampleType <- function(z, d, y, wc, nc0, nn0){
  if(z == 1 & d == 0){
    pc <- 0
  }else{
    if(z == 1 & d == 1){
      pc <- 1
    }else{
      if(y == 1){
        pc <- (wc * nc0) / (wc * nc0 +  (1 - wc) * nn0)
      }else{
        pc <- (wc * (1 - nc0)) / (wc * (1 - nc0) +  (1 - wc) * (1 - nn0))
      }
    }
  }
  rbinom(1, 1, prob = pc) # 1: Complier/ 0: Never taker
}
z = 0; d = 0; y = 0; wc = 0.8; nc0 = 0.9; nn0 = 0.05
SampleType(z = z, d = d, y = y, wc = wc, nc0 = nc0, nn0 = nn0)
Clat <- sapply(1:N, function(i){SampleType(z = Z[i], d = D[i], y = Y[i], wc = wc, nc0 = nc0, nn0 = nn0)})

# Gibbs sampler
a0 <- 1; b0 <- 1 # Hyperparameters beta priors
burnin <- 500; S <- 2000; tot <- S + burnin 
PosteriorDraws <- matrix(NA, tot, 5)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)

for(s in 1:tot){
  dataLat <- cbind(mydata, Clat)
  Nc011 <- sum(dataLat$Z == 0 & dataLat$Clat == 1 & dataLat$Y == 1)
  Nc010 <- sum(dataLat$Z == 0 & dataLat$Clat == 1 & dataLat$Y == 0)
  Nn001 <- sum(dataLat$Z == 0 & dataLat$Clat == 0 & dataLat$Y == 1)
  Nn000 <- sum(dataLat$Z == 0 & dataLat$Clat == 0 & dataLat$Y == 0)
  Nc <- sum(Clat == 1)
  Nn <- sum(Clat == 0)
  wc <- rbeta(1, 1 + Nc, 1 + Nn)
  nc1 <- rbeta(1, 1 + Nc111, 1 + Nc110)
  nn1 <- rbeta(1, 1 + Nn101, 1 + Nn100)
  nc0 <- rbeta(1, 1 + Nc011, 1 + Nc010)
  nn0 <- rbeta(1, 1 + Nn001, 1 + Nn000)
  Clat <- sapply(1:N, function(i){SampleType(z = Z[i], d = D[i], y = Y[i], wc = wc, nc0 = nc0, nn0 = nn0)})
  PosteriorDraws[s, ] <- c(wc, nc1, nc0, nn1, nn0)
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot)
LATE <- PosteriorDraws[keep, 2] - PosteriorDraws[keep, 3]
LATEmean <- mean(LATE)
LATEci <- quantile(LATE, c(0.025, 0.975))

# Imposing exclusion restrictions
PosteriorDrawsER <- matrix(NA, tot, 4)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)

for(s in 1:tot){
  dataLat <- cbind(mydata, Clat)
  Nc011 <- sum(dataLat$Z == 0 & dataLat$Clat == 1 & dataLat$Y == 1)
  Nc010 <- sum(dataLat$Z == 0 & dataLat$Clat == 1 & dataLat$Y == 0)
  Nn01 <- sum(dataLat$Clat == 0 & dataLat$Y == 1)
  Nn00 <- sum(dataLat$Clat == 0 & dataLat$Y == 0)
  Nc <- sum(Clat == 1)
  Nn <- sum(Clat == 0)
  wc <- rbeta(1, 1 + Nc, 1 + Nn)
  nc1 <- rbeta(1, 1 + Nc111, 1 + Nc110)
  nn <- rbeta(1, 1 + Nn01, 1 + Nn00)
  nc0 <- rbeta(1, 1 + Nc011, 1 + Nc010)
  Clat <- sapply(1:N, function(i){SampleType(z = Z[i], d = D[i], y = Y[i], wc = wc, nc0 = nc0, nn0 = nn)})
  PosteriorDrawsER[s, ] <- c(wc, nc1, nc0, nn)
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot)
LATEER <- PosteriorDrawsER[keep, 2] - PosteriorDrawsER[keep, 3]
LATEERmean <- mean(LATEER)
LATEERci <- quantile(LATEER, c(0.025, 0.975))


# Plot posterior distribution of ATE
hist(LATE, breaks = 40, freq = FALSE,
     main = "Posterior Distribution of CACE",
     xlab = "CACE", col = "lightblue", border = "white")
abline(v = LATEmean, col = "red", lwd = 2)
abline(v = LATEci, col = "darkgreen", lty = 2, lwd = 2)

legend("topright", legend = c("Posterior Mean", "95% Credible Interval"),
       col = c("red", "darkgreen"), lwd = 2, lty = c(1, 2),
       bty = "n", cex = 0.8)  # Smaller legend using cex

# Plot posterior distribution of ATE
hist(LATEER, breaks = 40, freq = FALSE,
     main = "Posterior Distribution of LATE: Exclusion restrictions",
     xlab = "LATE", col = "lightblue", border = "white")
abline(v = LATEERmean, col = "red", lwd = 2)
abline(v = LATEERci, col = "darkgreen", lty = 2, lwd = 2)

legend("topright", legend = c("Posterior Mean", "95% Credible Interval"),
       col = c("red", "darkgreen"), lwd = 2, lty = c(1, 2),
       bty = "n", cex = 0.8)  # Smaller legend using cex

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

#### Synthetic DiD data: parallel trends + no anticipation ####
rm(list = ls()); set.seed(10101)
library(ggplot2); library(dplyr); library(fastDummies)

# Parameters
N_per_group <- 200          # units per group
T_periods    <- 2           # keep 2x2 for clarity
tau_true     <- 1         # ATT
sigma_eps    <- 0.5         # noise SD

# Panel index
id  <- rep(1:(2*N_per_group), each = T_periods)
t   <- rep(1:T_periods, times = 2*N_per_group)

# Group: treated (D=1) vs control (D=0)
D   <- rep(c(rep(0, N_per_group), rep(1, N_per_group)), each = T_periods)

# Post indicator (t=2 is post)
post <- as.integer(t == 2)

# Unit fixed effects (random heterogeneity)
alpha_i <- rnorm(2*N_per_group, 0, 0.8)
alpha   <- alpha_i[id]

# Time effects (common shocks) – parallel trends built in here
phi_t <- c(0, -1.8)  # common decline from t=1 to t=2
phi   <- phi_t[t]

# No anticipation: effect is zero in pre, equals tau in post *only for treated*
treat_effect <- tau_true * (D * post)

# Outcome: Y_it(0) = alpha_i + phi_t + linear trend if you want; here just FE + shock
eps <- rnorm(length(id), 0, sigma_eps)
Y   <- alpha + phi + treat_effect + eps

did <- data.frame(id, t, D, post, Y)

# --- Quick check: cell means --------------------------------------------------
with(did, tapply(Y, list(Group = D, Time = t), mean))

# --- Plot group means over time (shows parallel trends + no anticipation) ----

plot_data <- did %>%
  group_by(D, t) %>%
  summarise(meanY = mean(Y), .groups = "drop") %>%
  mutate(Group = ifelse(D == 1, "Treated", "Control"))

ggp <- ggplot(plot_data, aes(x = t, y = meanY, group = Group, linetype = Group)) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_x_continuous(breaks = c(1, 2), labels = c("t = 1 (pre)", "t = 2 (post)")) +
  labs(x = "Time", y = "Mean outcome", title = "Synthetic DiD: Parallel Trends & No Anticipation") +
  theme_minimal(base_size = 12)
print(ggp)

# Bayesian inference: Model in differences
dY <- did$Y[did$t==2] - did$Y[did$t==1]
D  <- did$D[did$t==1]
post_fit <- MCMCpack::MCMCregress(dY ~ 1 + D, burnin = 100, mcmc = 1000)
tau_draws <- post_fit[, "D"]
quantile(tau_draws, c(.025,.5,.975))
quantile(post_fit[,1], c(.025,.5,.975))

# Bayesian inference: Model with interaction treatment x post period
post_fit1 <- MCMCpack::MCMCregress(Y ~ 1 + factor(id) + factor(t) + I(D * post), data = did, burnin = 100, mcmc = 1000)
tau_draws1 <- post_fit1[, "I(D * post)"]
quantile(tau_draws1, c(.025,.5,.975))

# Calculating the means at each level, and then ATT
# ---- Build 2x2 cell indicators (no intercept model) ----
did$c_pre  <- as.integer(did$D == 0 & did$t == 1)
did$c_post <- as.integer(did$D == 0 & did$t == 2)
did$t_pre  <- as.integer(did$D == 1 & did$t == 1)
did$t_post <- as.integer(did$D == 1 & did$t == 2)

# Sanity check: each row must belong to exactly one cell
stopifnot(all(did$c_pre + did$c_post + did$t_pre + did$t_post == 1))

# ---- Bayesian saturated cell-mean model (compatible priors) ----
# Y_it = mu_c,pre * 1{c,pre} + mu_c,post * 1{c,post} + mu_t,pre * 1{t,pre} + mu_t,post * 1{t,post} + eps
# No intercept: each coefficient IS a cell mean.
fit_cells <- MCMCpack::MCMCregress(
  Y ~ 0 + c_pre + c_post + t_pre + t_post,
  data  = did,
  burnin = 1000, mcmc = 10000, thin = 5
)

draws_cells <- as.matrix(fit_cells)
colnames(draws_cells)  # should be c("c_pre","c_post","t_pre","t_post", "sigma2")

# Posterior of the four cell means
mu_c_pre_draw  <- draws_cells[, "c_pre"]
mu_c_post_draw <- draws_cells[, "c_post"]
mu_t_pre_draw  <- draws_cells[, "t_pre"]
mu_t_post_draw <- draws_cells[, "t_post"]

# ATT as a linear combination of cell means:
# ATT = (mu_t,post - mu_t,pre) - (mu_c,post - mu_c,pre)
ATT_cells_draw <- (mu_t_post_draw - mu_t_pre_draw) - (mu_c_post_draw - mu_c_pre_draw)
quantile(ATT_cells_draw, c(.025, .5, .975))

# Frequentist
# --- Classic 2x2 DiD regression ----------------------------------------------
did$DID <- did$D * did$post
m_2x2 <- lm(Y ~ D + post + DID, data = did)
cat("\n2x2 DiD (ATE on treated) is coef on DID:\n")
print(summary(m_2x2)$coefficients["DID", , drop = FALSE])

# --- TWFE with unit and time FE + cluster-robust SEs (by unit) ---------------
# install.packages(c("fixest","sandwich","lmtest"), dependencies = TRUE)
library(fixest)    # fast FE estimation
library(sandwich)  # robust vcov
library(lmtest)    # coeftest

m_twfe <- feols(Y ~ D:post | id + t, data = did)  # same as Y ~ tau*(D*post) + unit FE + time FE
cat("\nTWFE estimate (coef on D:post):\n")
print(coeftest(m_twfe, vcov = vcovCL, cluster = ~ id))
cat("\nTrue ATT:", tau_true, "\n")

#### Comparing first difference specification with TWFE specification ####
rm(list = ls()); set.seed(10101)

sim_once <- function(N = 2000,
                     rhoX = 0.7,          # persistence of X
                     beta = 1.0,          # effect of X on Y
                     tau0 = 1.0,          # baseline ATT at X=0
                     tauX = 0.8,          # heterogeneity: ATT increases with X1
                     sel = 1.0,           # selection on X1 into treatment
                     sd_a = 1, sd_e = 1)  # unit FE and noise SDs
{
  #--- panel indices (2x2 DiD)
  id  <- 1:N
  t   <- rep(1:2, each = N)
  Post <- as.integer(t == 2)
  
  #--- pre/post covariate
  X1_pre  <- rnorm(N)
  X1_post <- rhoX * X1_pre + sqrt(1 - rhoX^2) * rnorm(N)
  
  #--- treatment depends on pre X
  D <- as.numeric(sel * X1_pre + rnorm(N) > 0)
  
  #--- heterog. treatment effect: tau_i = tau0 + tauX * X1_pre
  tau_i <- tau0 + tauX * X1_pre
  
  #--- unit FE and noise
  a  <- rnorm(N, 0, sd_a)
  e1 <- rnorm(N, 0, sd_e)
  e2 <- rnorm(N, 0, sd_e)
  
  #--- outcomes
  # levels model: Y_it(0) = a_i + phi_t + beta * X_it + eps_it
  phi <- c(0, -1)               # common trend (parallel under X)
  Y_pre  <- a + phi[1] + beta * X1_pre  + e1
  Y_post <- a + phi[2] + beta * X1_post + e2 + tau_i * D
  
  # truths to compare
  ATT_overall  <- mean(tau_i[D == 1])     # population ATT in this sample
  ATT_baseline <- tau0                    # effect when X1_pre = 0
  
  #--- stack panel
  Y   <- c(Y_pre, Y_post)
  X   <- c(X1_pre, X1_post)
  Di  <- rep(D, 2)
  idf <- factor(rep(id, 2))
  tt  <- factor(rep(1:2, each = N))
  dY  <- Y_post - Y_pre
  dX  <- X1_post - X1_pre
  
  #================= Estimation =================
  
  # (1) First-difference with ΔX: ΔY ~ D + ΔX
  fd <- lm(dY ~ D + dX)
  
  # (2) TWFE with interactions to identify baseline ATT:
  #     Y_it = α_i + φ_t + τ0*(D*Post) + γ*(Post*X1_pre) + δ*(D*Post*X1_pre) + ε_it
  # NOTE: D main effect is collinear with α_i and drops automatically.
  twfe <- lm(Y ~ idf + tt + I(D * (tt == 2)) +
               I((tt == 2) * X1_pre) + I(D * (tt == 2) * X1_pre))
  
  c(ATT_overall  = ATT_overall,
    ATT_baseline = ATT_baseline,
    FD_beta_D    = coef(fd)["D"],
    TWFE_tau0    = coef(twfe)["I(D * (tt == 2))"])
}

# Single run (illustrative)
res1 <- sim_once()
round(res1, 3)

# Monte Carlo to stabilize comparisons
R <- 100
res <- replicate(R, sim_once())
res <- t(res)

apply(res, 2, function(x) c(mean = mean(x), sd = sd(x))) |> round(3)


#### Working example: Staggered ATT ####
# Install packages if needed
library(did)
library(ggplot2)

#----------------------------
# 1. Simulate example data
#----------------------------
set.seed(123)

N <- 200   # number of individuals
T <- 5     # number of periods

# Balanced panel structure
df <- expand.grid(id = 1:N, year = 2000:(2000+T-1))

# Treatment adoption time (staggered: some treated earlier, some later, some never)
adopt_year <- sample(c(2002, 2003, 2004, 9999), N, replace = TRUE, prob = c(0.25, 0.25, 0.25, 0.25))
df$G <- adopt_year[df$id]  # group = first year treated (9999 for never treated)

# Outcome: baseline + trend + treatment effect if treated after adoption
tau <- 3
df$Y <- 10 + 0.5*(df$year - 2000) + rnorm(nrow(df))
df$Y <- df$Y + ifelse(df$year >= df$G & df$G != 9999, tau, 0)

#----------------------------
# 2. Estimate ATT using did
#----------------------------
out <- att_gt(
  yname = "Y",           # outcome
  tname = "year",        # time variable
  idname = "id",         # unit identifier
  gname = "G",           # first treatment period
  data = df,
  panel = TRUE,
  est_method = "dr"      # doubly robust (default)
)

summary(out)

#----------------------------
# 3. Event-study visualization
#----------------------------
es <- aggte(out, type = "dynamic")
summary(es)

ggdid(es) +
  labs(title = "Event-Study: ATT over time",
       x = "Time relative to treatment",
       y = "ATT")

#### Regression discontinuity design ####
rm(list = ls()); set.seed(10101)

# --- Simulation setup
N <- 3000
Zi <- sample(seq(-24, 24, 1), N, replace = TRUE)       # Running var
Z <- Zi - mean(Zi)
T <- as.integer(Z >= 0)                                # Policy indicator (0/1)
X <- sample(seq(85, 95, 1), N, replace = TRUE)         # Regressor

Wc <- cbind(1, Z, X)                                   # For compliers (includes Z)
W  <- cbind(1, X)                                      # For n/a (no Z effect)

B0c <- c(4.5, -0.2, 0.03)
B1c <- c(5.5, 0.4,  0.03)
B0n <- c(6.8,  -0.02)
B1a <- c(5.5,  -0.04)

s2c <- 0.1; s2n <- 0.15; s2a <- 0.2
Etac <- 0.7; Etan <- 0.15; Etaa <- 0.15
v <- 5

# --- Potential outcomes (Student-t noise scaled to match variances)
mu0c <- sqrt(s2c) * rt(N, v);  Y0c <- as.numeric(Wc %*% B0c + mu0c)
mu1c <- sqrt(s2c) * rt(N, v);  Y1c <- as.numeric(Wc %*% B1c + mu1c)
mu0n <- sqrt(s2n) * rt(N, v);  Y0n <- as.numeric(W  %*% B0n + mu0n)
mu1a <- sqrt(s2a) * rt(N, v);  Y1a <- as.numeric(W  %*% B1a + mu1a)

# --- Latent types: first n, then a, remaining are c
id    <- seq_len(N)
n_n   <- as.integer(round(Etan * N))
n_a   <- as.integer(round(Etaa * N))

id0n  <- sample(id, n_n, replace = FALSE)                  # never-takers
id_rem <- setdiff(id, id0n)
id1a  <- sample(id_rem, n_a, replace = FALSE)              # always-takers
idc   <- setdiff(id, c(id0n, id1a))                        # compliers

type <- rep(NA_character_, N)
type[id0n] <- "n"
type[id1a] <- "a"
type[idc]  <- "c"

# --- Realized treatment under monotonicity:
# n -> D=0, a -> D=1, c -> D=T
D <- integer(N)
D[type == "n"] <- 0
D[type == "a"] <- 1
D[type == "c"] <- T[type == "c"]

# --- 2x2 table of D vs T
tab_DT <- table(T = factor(T, levels = 0:1), D = factor(D, levels = 0:1))
print(tab_DT)

# --- Plot potential outcomes as functions of Z
# We'll show the simulated potential outcomes (not the realized Y).
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

df_plot <- data.frame(
  Z, X, T, 
  Y0c = Y0c, Y1c = Y1c, Y0n = Y0n, Y1a = Y1a
)

long <- df_plot |>
  select(Z, Y0c, Y1c, Y0n, Y1a) |>
  pivot_longer(cols = -Z, names_to = "outcome", values_to = "Y")

ggplot(long, aes(x = Z, y = Y)) +
  geom_point(alpha = 0.15, size = 1) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, span = 0.6) +
  facet_wrap(~ outcome, ncol = 2, scales = "free_y") +
  labs(title = "Potential outcomes as functions of Z",
       x = "Running variable Z", y = "Potential outcome") +
  theme_bw()

# Y
Y <- ifelse(type == "n", Y0n,
            ifelse(type == "a", Y1a,
                   ifelse(D == 0, Y0c, Y1c)))

# Plot Y vs Z
df_obs <- data.frame(Z = Z, Y = Y, D = factor(D), T = factor(T))

ggplot(df_obs, aes(x = Z, y = Y, color = D)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, span = 0.6) +
  labs(title = "Observed outcome Y vs running variable Z",
       x = "Running variable Z",
       y = "Observed outcome Y",
       color = "Treatment D") +
  theme_bw()

## Gibbs sampler ##
I00 <- which(T == 0 & D == 0)
I11 <- which(T == 1 & D == 1)
I01 <- which(T == 0 & D == 1)
I10 <- which(T == 1 & D == 0)
# Hyperparameters
a0c <- 1; a0n <- 1; a0a <- 1 # Hyperpar Dirichlet
v <- 5 # Student-t
a0 <- 0.01; d0 <- 0.01 # Inv-Gamma
# Prob complier | D=0
dens_t_locscale <- function(y, mu, sig2, v) {
  z <- (y - mu) / sqrt(sig2)
  dt(z, df = v) / sqrt(sig2)
}
# Prob complier | D=0
Pc0 <- function(b0c, sigma20c, etac, b0n, sigma20n, etan, i){
  p0c <- etac * dens_t_locscale(Y[i], mu = sum(Wc[i,]*b0c), sig2 = sigma20c, v = v)
  p0n <- etan * dens_t_locscale(Y[i], mu = sum(W[i, ]*b0n), sig2 = sigma20n, v = v)
  pc0 <- p0c / (p0c + p0n + 1e-12)
  return(pc0)
}
# Prob complier | D=1
Pc1 <- function(b1c, sigma21c, etac, b1a, sigma21a, etaa, i){
  p1c <- etac * dens_t_locscale(Y[i], mu = sum(Wc[i,]*b1c), sig2 = sigma21c, v = v)
  p1a <- etaa * dens_t_locscale(Y[i], mu = sum(W[i, ]*b1a), sig2 = sigma21a, v = v)
  pc1 <- p1c / (p1c + p1a + 1e-12)
  return(pc1)
}

# Function probability type
rtype <- function(alphas){
  types <- MCMCpack::rdirichlet(1, alphas)
  return(types)
}

PostLambda_i <- function(sig2, Beta, yi, hi, v){
  shape <- (v + 1)/2
  rate  <- (v + (yi - sum(hi * Beta))^2 / sig2)/2
  rgamma(1, shape = shape, rate = rate)
}

PostBeta <- function(sig2, lambda, y, H){
  k <- dim(H)[2]
  B0i <- solve(1000*diag(k))
  b0 <- rep(0, k)
  Bn <- solve(B0i + sig2^(-1)*t(H)%*%diag(lambda)%*%H)
  bn <- Bn%*%(B0i%*%b0 + sig2^(-1)*t(H)%*%diag(lambda)%*%y)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}

PostSig2 <- function(Beta, lambda, y, H){
  Ndk <- length(y)
  an <- a0 + Ndk 
  dn <- d0 + t(y - H%*%Beta)%*%diag(lambda)%*%(y - H%*%Beta)
  sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(sig2)
}

# MCMC parameter
burnin <- 1000; S <- 5000; tot <- S + burnin 
ETAs <- matrix(NA, tot, 3)
BETAS0C <- matrix(NA, tot, 3)
BETAS1C <- matrix(NA, tot, 3)
BETASA <- matrix(NA, tot, 2)
BETASN <- matrix(NA, tot, 2)
SIGMAS <- matrix(NA, tot, 4)
b0c <- rep(0, 3); b1c <- rep(0, 3); b1a <- rep(0, 2); b0n <- rep(0, 2)
sigma21c <- 0.1; sigma20c <- 0.1; sigma20n <- 0.1; sigma21a <- 0.1
EtacNew <- 0.5; EtanNew <- 0.25; EtaaNew <- 0.25
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  ProbC0 <- sapply(I00, function(i)Pc0(b0c, sigma20c, EtacNew, b0n, sigma20n, EtanNew, i = i))
  ProbC1 <- sapply(I11, function(i)Pc1(b1c, sigma21c, EtacNew, b1a, sigma21a, EtaaNew, i = i))
  Typec0 <- sapply(ProbC0, function(p) {sample(c("c", "n"), 1, prob = c(p, 1-p))})
  Typec1 <- sapply(ProbC1, function(p) {sample(c("c", "a"), 1, prob = c(p, 1-p))})
  typeNew <- rep(NA_character_, N)
  typeNew[I10] <- "n"
  typeNew[I01] <- "a"
  typeNew[I00]  <- Typec0
  typeNew[I11]  <- Typec1
  Nc <- sum(typeNew == "c")
  Na <- sum(typeNew == "a")
  Nn <- sum(typeNew == "n")
  anc <- a0c + Nc; ann <- a0n + Nn; ana <- a0a + Na
  alphas <- c(anc, ann, ana)
  EtasNew <- rtype(alphas = alphas)
  EtacNew <- EtasNew[1]; EtanNew <- EtasNew[2]; EtaaNew <- EtasNew[3]
  ETAs[s, ] <- EtasNew
  Lambda <- numeric(N)
  idx_n <- which(typeNew == "n")
  idx_a <- which(typeNew == "a")
  idx_c0 <- which(typeNew == "c" & D == 0)
  idx_c1 <- which(typeNew == "c" & D == 1)
  Lambda[idx_n]  <- sapply(idx_n,  function(i) PostLambda_i(sigma20n, b0n, Y[i],  W[i, ], v))
  Lambda[idx_a]  <- sapply(idx_a,  function(i) PostLambda_i(sigma21a, b1a, Y[i],  W[i, ], v))
  Lambda[idx_c0] <- sapply(idx_c0, function(i) PostLambda_i(sigma20c, b0c, Y[i], Wc[i, ], v))
  Lambda[idx_c1] <- sapply(idx_c1, function(i) PostLambda_i(sigma21c, b1c, Y[i], Wc[i, ], v))
  b1a <- PostBeta(sig2 = sigma21a, lambda = Lambda[idx_a], y = Y[idx_a], H = W[idx_a,])
  b0n <- PostBeta(sig2 = sigma20n, lambda = Lambda[idx_n], y = Y[idx_n], H = W[idx_n,])
  b0c <- PostBeta(sig2 = sigma20c, lambda = Lambda[idx_c0], y = Y[idx_c0], H = Wc[idx_c0,])
  b1c <- PostBeta(sig2 = sigma21c, lambda = Lambda[idx_c1], y = Y[idx_c1], H = Wc[idx_c1,])
  BETASA[s, ] <- b1a
  BETASN[s, ] <- b0n
  BETAS0C[s, ] <- b0c
  BETAS1C[s, ] <- b1c 
  sigma21a <- PostSig2(Beta = b1a, lambda = Lambda[idx_a], y = Y[idx_a], H = W[idx_a,])
  sigma20n <- PostSig2(Beta = b0n, lambda = Lambda[idx_n], y = Y[idx_n], H = W[idx_n,])
  sigma20c <- PostSig2(Beta = b0c, lambda = Lambda[idx_c0], y = Y[idx_c0], H = Wc[idx_c0,])
  sigma21c <- PostSig2(Beta = b1c, lambda = Lambda[idx_c1], y = Y[idx_c1], H = Wc[idx_c1,])
  SIGMAS[s, ] <- c(sigma21a, sigma20n, sigma20c, sigma21c)
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0),"% done"))
}
close(pb)
keep <- seq((burnin+1), tot)
summary(coda::mcmc(BETASA[keep, ]))
summary(coda::mcmc(BETASN[keep, ]))
summary(coda::mcmc(BETAS0C[keep, ]))
summary(coda::mcmc(BETAS1C[keep, ]))
summary(coda::mcmc(SIGMAS[keep, ]))
LATE <- coda::mcmc(BETAS1C[keep, 1] - BETAS0C[keep, 1])
summary(LATE)
plot(LATE)
# Extract samples as numeric
late_draws <- as.numeric(LATE)
# Posterior mean and 95% CI
LATEmean <- mean(late_draws)
LATEci   <- quantile(late_draws, c(0.025, 0.975))
# Plot posterior density
df <- data.frame(LATE = late_draws)

ggplot(df, aes(x = LATE)) +
  geom_density(fill = "skyblue", alpha = 0.5, color = "blue") +
  geom_vline(xintercept = LATEmean, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = LATEci, color = "black", linetype = "dotted", linewidth = 1) +
  labs(title = "Posterior distribution of LATE",
       subtitle = paste0("Mean = ", round(LATEmean,3),
                         " | 95% CI = [", round(LATEci[1],3), ", ", round(LATEci[2],3), "]"),
       x = "LATE", y = "Density") +
  theme_bw(base_size = 14)

#### Sample selection models #####
rm(list = ls()); set.seed(10101)
library(doParallel); library(snow)
N <- 1000
w1 <- rbinom(N, 1, 0.5) # rnorm(N, 0, sigExo); 
delta <- c(1, 1, -1); beta <- c(2,1,1)
zx <- MASS::mvrnorm(N, mu = rep(0, 2), matrix(c(1,0.7,0.7,1),2,2))
z1 <- zx[,1]
Z <- cbind(1, z1, w1)
sig12 <- 0.8
sig11 <- 1.2
SIGMA <- matrix(c(sig11, sig12, sig12, 1), 2, 2)
E <- MASS::mvrnorm(N, mu = rep(0, 2), SIGMA)
cl <- Z%*%delta + E[,2]
c <- cl > 0
table(c)
x1 <- zx[,2]; X <- cbind(1, x1, w1)
y <- X%*%beta + E[,1]
y[c==0] <- NA
# Hyperparameters
b0 <- rep(0, 6); B0 <- 1000*diag(6); B0i <- solve(B0)
a0 <- 0.001; d0 <- 0.001
s0 <- 0; S0 <- 1000; S0i <- 1/S0
# Location
idc1 <- which(c==1)
nc1 <- length(idc1)
PostThetaNew <- function(Sigma, clat){
  J <- matrix(c(0,0,0,1),2,2)
  WW <- matrix(0, 6, 6)
  Wy <- matrix(0, 6, 1)
  for(i in 1:N){
    if(i %in% idc1){
      yclat <- c(y[i], clat[i])
      Auxi <- solve(Sigma)
    }else{
      yclat <- c(0, clat[i])
      Auxi <- J
    }
    Wi <- as.matrix(Matrix::bdiag(X[i,], Z[i,]))
    WWi <- Wi%*%Auxi%*%t(Wi)
    Wyi <- Wi%*%Auxi%*%yclat
    WW <- WW + WWi
    Wy <- Wy + Wyi
  }
  Bn <- solve(B0i + WW)
  bn <- Bn%*%(B0i%*%b0 + Wy)
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}
# clat <- cl; Sigma <- SIGMA[1:2,1:2]
# PostThetaNew(Sigma = Sigma, clat = clat)

PostOmega11 <- function(theta, sig12, clat){
  # sig12 <- SIGMA[1,2]; theta <- c(beta, delta); clat <- cl
  an <- a0 + nc1
  mui <- y[idc1] -X[idc1, ]%*%theta[1:3] - sig12*(clat[idc1] - Z[idc1,]%*%theta[4:6])
  dn <- d0 + t(mui)%*%mui
  omega11 <- LaplacesDemon::rinvgamma(1, an/2, dn/2)
  return(omega11)
}
# sig12 <- SIGMA[1,2]; theta <- c(beta, delta); clat <- cl
# PostOmega11(theta, sig12, clat)

PostSig12 <- function(omega11, theta, clat){
  Sn <- (omega11^(-1)*sum((clat[idc1] - Z[idc1,]%*%theta[4:6])^2) + S0i)^(-1)
  sn <- Sn*(omega11^(-1)*sum((clat[idc1] - Z[idc1,]%*%theta[4:6])*(y[idc1] - X[idc1,]%*%theta[1:3])) + s0*S0i)
  sig12 <- rnorm(1, sn, sd = Sn^0.5)
  return(sig12)
}
# omega11 <- SIGMA[1,1] - SIGMA[1,2]^2; theta <- c(beta, delta); clat <- cl
# PostSig12(omega11, theta, cl)
# PostOmega11(theta, sig12, clat)+PostSig12(omega11, theta, cl)^2

PostClat <- function(theta, sig12, omega11, i){
  if(i %in% idc1){
    mu <- Z[i,]%*%theta[4:6] + (sig12/(omega11+sig12^2))*(y[i] - X[i,]%*%theta[1:3])
    sig2 <- omega11/(omega11+sig12^2)
    clat <- EnvStats::rnormTrunc(1, mean = mu, sd = sig2^0.5, min = 0, max = Inf)
  }else{
    mu <- Z[i,]%*%theta[4:6]
    clat <- EnvStats::rnormTrunc(1, mean = mu, sd = 1, min = -Inf, max = 0)
  }
  return(clat)
}
# i <- 3
# omega11 <- SIGMA[1,1] - SIGMA[1,2]^2; sig12 <- SIGMA[1,2]; theta <- c(beta, delta)
# PostClat(theta, sig12, omega11, i)

# Sampler
S <- 1500
burnin <- 500
thin <- 2
keep <- seq(burnin+thin, S, thin)
PostThetasDraws <- matrix(NA, S, 6)
PostSigmaDraws <- matrix(NA, S, 2)

# Initial conditions
Thetap <- rep(0, 6)
Sigmap <- diag(2) 
Sig12p <- 0
Omega11p <- 1
for(i in 1:N){
  if(c[i] == 0){
    LatPost <- EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = -Inf, max = 0)
  }else{
    LatPost <- EnvStats::rnormTrunc(1, mean = 0, sd = 1, min = 0, max = Inf)
  }
}

#### Parallel code ####
cn <- detectCores() 
ClusterHope <- makeCluster(cn, type = "SOCK")
registerDoParallel(ClusterHope)
clusterExport(ClusterHope, list("Z", "X", "c", "y", "N", "idc1", "nc1", 
                                "PostClat","Thetap", "Sig12p", "Omega11p"))
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
for(rep in 1:S){
  LatPost <- t(parSapply(ClusterHope, 1:N, function(i){PostClat(theta = Thetap, sig12 = Sig12p, omega11 = Omega11p, i)}))
  Thetap <- PostThetaNew(Sigma = Sigmap, clat = LatPost)
  Omega11p <- PostOmega11(theta = Thetap, sig12 = Sig12p, clat = LatPost)
  Sig12p <- PostSig12(omega11 = Omega11p, theta = Thetap, clat = LatPost)
  Sigmap <- matrix(c(Omega11p+Sig12p^2, Sig12p, Sig12p, 1),2,2)
  PostThetasDraws[rep,] <- Thetap
  PostSigmaDraws[rep, ] <- c(Omega11p+Sig12p^2, Sig12p)
  clusterExport(ClusterHope, list("Thetap", "Sig12p", "Omega11p"))
  setWinProgressBar(pb, rep, title=paste( round(rep/S*100, 0),"% done"))
}
stopCluster(ClusterHope)
close(pb)
thetaHat <- coda::mcmc(PostThetasDraws[keep,])
PostDrawsSigma <- coda::mcmc(PostSigmaDraws[keep,])
summary(thetaHat)
summary(PostDrawsSigma)

RegNOsel <- MCMCpack::MCMCregress(y ~ X - 1)
summary(RegNOsel)

# --- Inputs you already have:
# RegNOsel <- MCMCpack::MCMCregress(y ~ X - 1)
# thetaHat <- coda::mcmc(PostThetasDraws[keep, ])   # takes selection into account

library(coda)
library(ggplot2)

# 1) Extract the 2nd coefficient draws from each object
# For MCMCpack::MCMCregress, coefficients are first columns of the mcmc matrix.
# We'll grab column 2 explicitly. Adjust if your 2nd coefficient has a name you prefer.
beta2_nosel <- as.numeric(RegNOsel[, 2])

# For thetaHat (coda mcmc), also grab the 2nd parameter.
beta2_sel   <- as.numeric(thetaHat[, 2])

# 2) Put into a long data frame for ggplot
df <- rbind(
  data.frame(value = beta2_nosel, model = "No selection"),
  data.frame(value = beta2_sel,   model = "With selection")
)

# 3) Compute quick summaries
summ <- function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975)))
cat("\nPosterior summaries (2nd coefficient)\n")
cat("No selection : ", paste(round(summ(beta2_nosel), 4), collapse = "  "), "\n")
cat("With selection: ", paste(round(summ(beta2_sel),   4), collapse = "  "), "\n\n")

# 4) Plot both posteriors + true value line at 1
ggplot(df, aes(x = value, fill = model, color = model)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = 2) +
  labs(x = "Coefficient (2nd parameter)", y = "Posterior density",
       title = "Posterior of 2nd Coeff: No-Selection vs Selection Models",
       subtitle = "Dashed line = population value (1)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

#### BETEL ####

# To install betel package
# Go to: https://apps.olin.wustl.edu/faculty/chib/rpackages/betel/
# Download: betel_1.0.zip (windows) or betel_1.0.tgz (mac)
# Then Rstudio -> Tools -> Install Package -> and changing the Install from option to
# Package Archive file. Then, scroll to the location where the downloaded package was
# saved (typically the download files folder), and select the file.

#### BETEL: Measurement Error #####
rm(list = ls()); set.seed(10101)
library(betel); library(ucminf)
# Simulate data
N <- 2000; d <- 2; k <- 2
gamma <- 0.5; beta <- c(1, 1.2)
# Mixture
mum1 <- 1/2; mum2 <- -1/2
mu1 <- rnorm(N, mum1, 0.5); mu2 <- rnorm(N, mum2, 1.2)
mu <- sapply(1:N, function(i){sample(c(mu1[i], mu2[i]), 1, prob = c(0.5, 0.5))})
e <- rnorm(N)
z <- rnorm(N) # Instrument
xlat <- gamma*z + e # Unobserved regressor
nu <- rnorm(N) # Measurement error
x <- xlat + nu # Observed regressor
Xlat <- cbind(1, xlat)
y <- Xlat%*%beta + mu
dat <- cbind(1, x, z) # Data
# Function g_i by row in BETEL
gfunc <- function(psi = psi, y = y, dat = dat) {
  X <- dat[,1:2]
  e <- y - X %*% psi
  E <- e %*% rep(1,d)
  Z <- dat[,c(1,3)]
  G <- E * Z;
  return(G)
}
nt <- round(N * 0.1, 0); # training sample size for prior
psi0 <- lm(y[1:nt]~x[1:nt])$coefficients # Starting value of psi = (theta, v), v is the slack parameter in CSS (2018)
names(psi0) <- c("alpha","beta")
psi0_ <- as.matrix(psi0) # Prior mean of psi 
Psi0_ <- 5*rep(1,k) # Prior dispersions of psi
lam0 <- .5*rnorm(d) # Starting value of lambda
nu <- 2.5 # df of the prior student-t
nuprop <- 15 # df of the student-t proposal
n0 <- 1000 # burn-in
m <- 10000 # iterations beyond burn-in
# MCMC ESTIMATION BY THE CSS (2018) method
psim <- betel::bayesetel(gfunc = gfunc, y = y[-(1:nt)], dat = dat[-(1:nt),], psi0 = psi0, lam0 = lam0, psi0_ = psi0_, Psi0_ = Psi0_, nu = nu, nuprop = nuprop,
                         controlpsi = list(maxiterpsi = 50, mingrpsi = 1.0e-8), #  list of parameters in maximizing likelihood over psi
                         controllam = list(maxiterlam = 50, # list of parameters in minimizing dual over lambda
                                           mingrlam = 1.0e-7),
                         n0 = n0, m = m)

MCMCexg <- MCMCpack::MCMCregress(y ~ x, burnin = n0, mcmc = m)
Data <- list(y = c(y), x = x, z = matrix(z, N, 1), w = matrix(rep(1, N), N, 1))
Mcmc <- list(R = m)
MCMCivr <- bayesm::rivGibbs(Data, Mcmc = Mcmc)

dfplot <- data.frame(betel = psim[,2], iv = MCMCivr[["betadraw"]], exo = MCMCexg[,2])
colnames(dfplot) <- c("betel", "iv", "exo")

library(tidyr)
library(dplyr)
library(ggplot2)

df_long <- dfplot |>
  pivot_longer(everything(), names_to = "Method", values_to = "Posterior") |>
  mutate(Method = factor(Method,
                         levels = c("betel","iv","exo"),
                         labels = c("betel","rivGibbs","MCMCregress")))

ggplot(df_long, aes(x = Posterior, color = Method, fill = Method)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  geom_vline(xintercept = 1.2, linetype = "dashed", linewidth = 1, color = "black") +
  labs(
    title = "Posterior Densities with Population Value",
    x = expression(beta[1]),
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

summary(AER::ivreg(y ~ x | z))
summary(lm(y ~ x))

############# Instrumental Variable Quantile Regression: Simulation ##############
rm(list = ls()); set.seed(10101)
library(quantreg)
### Siulation ####
# Asymmetric Laplace for general τ in (0,1) with location 0, scale > 0
# Normal–Exponential mixture (Kozumi–Kobayashi, 2011)
rALD <- function(n, tau) {
  stopifnot(tau > 0, tau < 1)
  theta <- (1 - 2*tau) / (tau * (1 - tau))
  psi   <- 2 / (tau * (1 - tau))
  v  <- rexp(n, rate = 1)
  z  <- rnorm(n)
  mu <- theta * v + sqrt(psi * v) * z
  return(mu)
}

## --- Simulator for a given τ ---
simulate_qr <- function(n = 5000, tau = 0.5,
                        beta0 = 0.7,
                        beta1_tau = 1.5,
                        beta2_tau = 1,
                        gamma = 1,
                        delta = 0.5) {
  # Nonnegative regressor stabilizes intercept & preserves τ-ordering of slopes
  x <- runif(n, 0, 1)
  # optional: anchor a portion at 0 to pin intercept
  x[sample.int(n, size = round(0.1*n))] <- 0
  mu <- rALD(n, tau = tau) # ALD error
  z <- rnorm(n) # Instrument
  d <- gamma * z + delta * mu
  y <- beta0 + beta1_tau * x + beta2_tau * d + mu
  df <- data.frame(y, x, d, z)
}

tau <- 0.75
df <- simulate_qr(tau = tau, beta0 = 1, beta1_tau = 1,
                  beta2_tau = 1, gamma = 1, delta = 1, n = 5000)
Reg <- MCMCpack::MCMCquantreg(y~x+d, data = df, tau = tau)
summary(Reg)
summary(rq(y ~ x + d, tau = tau, data = df))

LossFunct <- function(par, tau, y, z, x, d){
  n <- length(y)
  X <- cbind(1, x, d)
  Z <- cbind(1, x, z)
  Ind <- as.numeric(y <= X%*%par) 
  gn <- colMeans((tau - Ind) * Z)
  Wni <- lapply(1:n, function(i) {Z[i,] %*% t(Z[i,])}) 
  Wn <- 1 / (tau * (1 - tau)) * solve(Reduce("+", Wni)/n)
  Ln <- - 0.5 * n * t(gn) %*% Wn %*% gn
  return(Ln)
}

par0 <- rnorm(3)
LossFunct(par = par0, tau = tau, y = df$y, z = df$z, x = df$x, d = df$d)

# ----- MH using Ln -----
k <- 3
b0 <- rep(0, k); B0 <- 1000*diag(k)
S <- 5000; burnin <- 6000; thin <- 5; tot <- S + burnin
BETA <- matrix(NA, tot, k); accept <- logical(tot)
step <- 0.05

BETA[1,] <- par0
LL <- LossFunct(par = BETA[1,], tau = tau, y = df$y, z = df$z, x = df$x, d = df$d)

pb <- txtProgressBar(min=0, max=tot, style=3)
for(s in 2:tot){
  cand <- BETA[s-1,] + rnorm(k, 0, step)
  LLc  <- LossFunct(par = cand, tau = tau, y = df$y, z = df$z, x = df$x, d = df$d)
  priorRat <- mvtnorm::dmvnorm(cand, b0, B0, log=TRUE) -
    mvtnorm::dmvnorm(BETA[s-1,], b0, B0, log=TRUE)
  loga <- (LLc - LL) + priorRat
  if (is.finite(loga) && log(runif(1)) <= loga) {
    BETA[s,] <- cand; LL <- LLc; accept[s] <- TRUE
  } else {
    BETA[s,] <- BETA[s-1,]; accept[s] <- FALSE
  }
  if (s <= burnin && s %% 200 == 0) {          # gentle adaptation
    acc <- mean(accept[(s-199):s])
    if (acc > 0.4) step <- step * 1.25
    if (acc < 0.15) step <- step / 1.25
  }
  setTxtProgressBar(pb, s)
}
close(pb)

cat("\nAcceptance rate:", mean(accept), "\n")

keep <- seq(burnin, tot, by = thin)
post <- BETA[keep, , drop=FALSE]   # posterior draws in scaled space

colnames(post) <- c("beta0","beta_x","beta_d")

# Posterior summaries
post_mean <- colMeans(post)
post_ci   <- apply(post, 2, quantile, c(0.025, 0.975))

round(post_mean, 3)
round(post_ci, 3)

plot(coda::mcmc(post))

