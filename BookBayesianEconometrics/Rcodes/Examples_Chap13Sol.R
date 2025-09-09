#### DAG: CIA  ####
library(dagitty)
library(ggdag)

Gd = dagitty('dag{
X [pos="-1,1"]
D [exposure, pos="0,0"]
Y [outcome, pos="1,1"]
X -> D
X -> Y
D -> Y
}')

ggdag(Gd) +  theme_dag()

isAcyclic(Gd)

for( n in names(Gd) ){
  for( m in children(Gd,n) ){
    a <- adjustmentSets(Gd, n, m, effect = c("direct"), type = c("minimal"))
    if( length(a) > 0 ){
      cat("The effect ",n,"->",m,
          " is identifiable by controlling for:\n",sep="")
      print( a, prefix=" * " )
    }
  }
}

#### DAG: Collider bias  ####
library(dagitty)
library(ggdag)

Gd = dagitty('dag{
X [pos="-1,1"]
D [exposure, pos="0,0"]
Y [outcome, pos="1,1"]
C [pos="0,0.5"]
X -> D
X -> Y
X -> C
D -> Y
D -> C
}')

ggdag(Gd) +  theme_dag()

isAcyclic(Gd)

adjustmentSets(Gd, exposure = "D", outcome = "Y")

for( n in names(Gd) ){
  for( m in children(Gd,n) ){
    a <- adjustmentSets(Gd, n, m, effect = c("direct"), type = c("minimal"))
    if( length(a) > 0 ){
      cat("The effect ",n,"->",m,
          " is identifiable by controlling for:\n",sep="")
      print( a, prefix=" * " )
    }
  }
}

#### DAG: Instruments  ####
library(dagitty)
library(ggdag)

Gd = dagitty('dag{
U [latent, pos="-1,1"]
D [exposure, pos="0,0"]
Y [outcome, pos="1,1"]
Z [pos="0,0.5"]
U -> D
U -> Y
D -> Y
Z -> D
}')

ggdag(Gd) +  theme_dag()

isAcyclic(Gd)

instrumentalVariables(Gd)

adjustmentSets(Gd, exposure = "D", outcome = "Y")

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

######## Instrumental variables: 401k participation #############
rm(list = ls()); set.seed(10101)
library(ggplot2); library(bayesm)

# Load data
mydata <- read.csv("https://raw.githubusercontent.com/BEsmarter-consultancy/BSTApp/refs/heads/master/DataApp/401k.csv",
                   sep = ",", header = TRUE, quote = "")
# Attach variables
attach(mydata)
y <- net_tfa/1000                  # Outcome: net financial assets
x <- as.vector(p401)          # Endogenous regressor: participation
w <- as.matrix(cbind(age, inc, fsize, educ, marr, twoearn, db, pira, hown))  # Exogenous regressors with intercept
z <- as.matrix(e401)          # Instrument: eligibility (NO intercept here)

# specify data input and mcmc parameters
Data = list(); 
Data$z <- z
Data$x <- x
Data$y <- y
Data$w <- w

Mcmc = list()
Mcmc$maxuniq <- 100
Mcmc$R <- 5000
keep <- seq((1000+1), Mcmc$R)

out <- rivDP(Data=Data, Mcmc=Mcmc)

# Convert to data frame for ggplot
df_ITT <- data.frame(ITT = out[["betadraw"]][kepp])

# Plot posterior distribution of treatment effect
ggplot(df_ITT, aes(x = ITT)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = mean(ITT), color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Posterior Distribution of 401(k) ITT: Dirichlet process",
    x = expression(beta["ITT"]),
    y = "Density"
  ) +
  theme_minimal(base_size = 14)

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

#### BETEL ####

# To install betel package
# Go to: https://apps.olin.wustl.edu/faculty/chib/rpackages/betel/
# Download: betel_1.0.zip (windows) or betel_1.0.tgz (mac)
# Then Rstudio -> Tools -> Install Package -> and changing the Install from option to
# Package Archive file. Then, scroll to the location where the downloaded package was
# saved (typically the download files folder), and select the file.

#### BETEL: Omission relevant correlated regressor #####
rm(list = ls()); set.seed(10101)
library(betel); library(ucminf)
# Simulate data
N <- 2000; d <- 2; k <- 2
gamma <- 0.5; beta <- c(1, 1.2, -0.7); rho <- 0.7
# Mixture
mum1 <- 1/2; mum2 <- -1/2
mu1 <- rnorm(N, mum1, 0.5); mu2 <- rnorm(N, mum2, 1.2)
mu <- sapply(1:N, function(i){sample(c(mu1[i], mu2[i]), 1, prob = c(0.5, 0.5))})
z <- rnorm(N) # Instrument
E <- MASS::mvrnorm(n = N, mu = c(0, 0), Sigma = matrix(c(1,rho,rho,1),2,2))
x <- gamma*z + E[,1] # Observed regressor
w <- rnorm(N) + E[,2] # Unobserved correlated regressor
X <- cbind(1, x, w)
y <- X%*%beta + mu
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
psim2 <- betel::bayesetel(gfunc = gfunc, y = y[-(1:nt)], dat = dat[-(1:nt),], psi0 = psi0, lam0 = lam0, psi0_ = psi0_, Psi0_ = Psi0_, nu = nu, nuprop = nuprop,
                         controlpsi = list(maxiterpsi = 50, mingrpsi = 1.0e-8), #  list of parameters in maximizing likelihood over psi
                         controllam = list(maxiterlam = 50, # list of parameters in minimizing dual over lambda
                                           mingrlam = 1.0e-7),
                         n0 = n0, m = m)

MCMCreg2 <- MCMCpack::MCMCregress(y~x)
dfplot <- data.frame(betel1 = psim1[,2], betel2 = psim2[,2], Exo = MCMCreg1[,2], Exo1 = MCMCreg2[,2])
colnames(dfplot) <- c("betel0.0", "betel0.7", "MCMCreg0.0", "MCMCreg0.7")
library(tidyr)
library(dplyr)
library(ggplot2)

df_long <- dfplot |>
  pivot_longer(everything(), names_to = "Method", values_to = "Posterior") |>
  mutate(Method = factor(Method,
                         levels = c("betel0.0","betel0.7", "MCMCreg0.0", "MCMCreg0.7"),
                         labels = c("BETEL: 0","BETEL: 0.7","Exo: 0", "Exo: 0.7")))

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

##################### BETEL: Demand vs Supply ######################
# Simulation
rm(list = ls()); set.seed(12345)
# Population parameters demand
B1 <- 5; B2 <- -0.5; B3 <- 0.8; B4 <- -0.4; B5 <- 0.7; SD <- 0.5
# Population parameters supply
A1 <- -2; A2 <- 0.5; A3 <- -0.4; SS <- 0.5
# Reduced form parameters
P0 <- (A1-B1)/(B2-A2); P2 <- -B3/(B2-A2); P3 <- -B4/(B2-A2); P1 <- A3/(B2-A2); P4 <- -B5/(B2-A2)
T0 <- B1+B2*P0; T2 <- B3+B2*P2; T3 <- B4+B2*P3; T1 <- B2*P1; T4 <- B5+B2*P4;
n <- 5000
ED <- rnorm(n, 0, SD); ES <- rnorm(n, 0, SS)
VP <- (ES-ED)/(B2-A2); UQ <- B2*VP+ED
y <- rnorm(n, 10, 1); pc <- rnorm(n, 5, 1); er <- rnorm(n, 15, 1); ps <- rnorm(n, 5, 1);
p <- P0+P1*er+P2*y+P3*pc+P4*ps+VP
q <- T0+T1*er+T2*y+T3*pc+T4*ps+UQ

# Supply
p_s  <- as.numeric(scale(p,  TRUE, TRUE))
er_s <- as.numeric(scale(er, TRUE, TRUE))

X <- cbind(1, p_s, er_s)  # k = 3  (STRUCTURAL: intercept, p, er)
Z <- cbind(1,
           as.numeric(scale(y,  TRUE, TRUE)),
           as.numeric(scale(pc, TRUE, TRUE)),
           as.numeric(scale(ps, TRUE, TRUE)),
           er_s)          # EXCLUDED instruments y,pc,ps, plus er

k <- ncol(X); d <- ncol(Z)

# ----- objective (stable) -----
lambdafunc <- function(lambda, theta, q, X, Z){
  e <- as.numeric(q - X %*% theta); if (!all(is.finite(e))) return(1e300)
  G <- Z * e
  eta <- as.numeric(G %*% lambda); if (!all(is.finite(eta))) return(1e300)
  m <- max(eta); m + log(mean(exp(eta - m)))
}

# ----- Profile ETEL log-likelihood (use this in MH) -----
loglik_etel <- function(theta, q, X, Z){
  e <- as.numeric(q - X %*% theta); if (!all(is.finite(e))) return(-Inf)
  G <- Z * e
  op <- optim(par = rep(0, ncol(G)), fn = function(l){
    eta <- as.numeric(G %*% l); m <- max(eta)
    m + log(mean(exp(eta - m)))
  }, method="BFGS", control=list(maxit=2000))
  if (!is.finite(op$value)) return(-Inf)
  lam <- op$par
  eta <- as.numeric(G %*% lam)
  sum(eta) - length(e) * log(sum(exp(eta)))   # ∑ log p_i(θ)
}

theta0 <- coef(lm(q ~ p_s + er_s))
# ----- MH using profile ETEL -----
b0 <- rep(0, k); B0 <- 1000*diag(k)
S <- 10000; burnin <- 2000; thin <- 5; tot <- S + burnin
BETA <- matrix(NA, tot, k); accept <- logical(tot)
step <- c(0.05, 0.02, 0.02)                  # preconditioned steps

BETA[1,] <- theta0
LL <- loglik_etel(BETA[1,], q, X, Z)

pb <- txtProgressBar(min=0, max=tot, style=3)
for(s in 2:tot){
  cand <- BETA[s-1,] + rnorm(k, 0, step)
  LLc  <- loglik_etel(cand, q, X, Z)
  priorRat <- mvtnorm::dmvnorm(cand, b0, B0, log=TRUE) -
    mvtnorm::dmvnorm(BETA[s-1,], b0, B0, log=TRUE)
  loga <- (LLc - LL) + priorRat
  if (is.finite(loga) && log(runif(1)) <= loga) {
    BETA[s,] <- cand; LL <- LLc; accept[s] <- TRUE
  } else {
    BETA[s,] <- BETA[s-1,]; accept[s] <- FALSE
  }
  if (s %% 200 == 0) {          # gentle adaptation
    acc <- mean(accept[(s-199):s])
    if (acc > 0.4) step <- step * 1.25
    if (acc < 0.15) step <- step / 1.25
  }
  setTxtProgressBar(pb, s)
}
close(pb)

cat("\nAcceptance rate:", mean(accept), "\n")

keep <- seq(burnin, tot, by = thin)
post_s <- BETA[keep, , drop=FALSE]   # posterior draws in scaled space

## ---------- Back-transform each draw to original scale ----------
# If q = β0 + βp * (p - μp)/σp + βer * (er - μer)/σer + u
# Then:
# αp   = βp / σp
# αer  = βer / σer
# α0   = β0 - βp * μp/σp - βer * μer/σer
p_mu <- mean(p);   p_sd <- sd(p)
er_mu <- mean(er); er_sd <- sd(er)
alpha_draws <- cbind(
  alpha0  = post_s[,1] - post_s[,2] * (p_mu/p_sd) - post_s[,3] * (er_mu/er_sd),
  alphap  = post_s[,2] / p_sd,
  alphaer = post_s[,3] / er_sd
)

## ---------- Summaries ----------
summ <- function(x) c(mean=mean(x), sd=sd(x), quantile(x, c(.025,.25,.5,.75,.975)))
cat("\nPosterior (original scale):\n")
print(rbind(
  alpha0  = summ(alpha_draws[,"alpha0"]),
  alphap  = summ(alpha_draws[,"alphap"]),
  alphaer = summ(alpha_draws[,"alphaer"])
))
cat("\nTrue (original scale):\n")
print(c(alpha0=A1, alphap=A2, alphaer=A3))
ElastSupPrice <- alpha_draws[,2]

# Demand
y_s <- as.numeric(scale(y))
pc_s <- as.numeric(scale(pc))
ps_s <- as.numeric(scale(ps))
Xdem <- cbind(1, p_s, y_s, pc_s, ps_s) 
k <- ncol(Xdem); d <- ncol(Z)

theta0 <- coef(lm(q ~ p_s + y_s + pc_s + ps_s))
# ----- MH over θ using profile ETEL -----
b0 <- rep(0, k); B0 <- 1000*diag(k)
S <- 10000; burnin <- 2000; thin <- 5; tot <- S + burnin
BETAdem <- matrix(NA, tot, k); accept <- logical(tot)
step <- c(0.05, 0.02, 0.02)                  # preconditioned steps

BETAdem[1,] <- theta0
LL <- loglik_etel(BETAdem[1,], q, X = Xdem, Z)

pb <- txtProgressBar(min=0, max=tot, style=3)
for(s in 2:tot){
  cand <- BETAdem[s-1,] + rnorm(k, 0, step)
  LLc  <- loglik_etel(cand, q, X = Xdem, Z)
  priorRat <- mvtnorm::dmvnorm(cand, b0, B0, log=TRUE) -
    mvtnorm::dmvnorm(BETAdem[s-1,], b0, B0, log=TRUE)
  loga <- (LLc - LL) + priorRat
  if (is.finite(loga) && log(runif(1)) <= loga) {
    BETAdem[s,] <- cand; LL <- LLc; accept[s] <- TRUE
  } else {
    BETAdem[s,] <- BETAdem[s-1,]; accept[s] <- FALSE
  }
  if (s %% 200 == 0) {          # gentle adaptation
    acc <- mean(accept[(s-199):s])
    if (acc > 0.4) step <- step * 1.25
    if (acc < 0.15) step <- step / 1.25
  }
  setTxtProgressBar(pb, s)
}
close(pb)

cat("\nAcceptance rate:", mean(accept), "\n")

keep <- seq(burnin, tot, by = thin)
postDem_s <- BETAdem[keep, , drop=FALSE]   # posterior draws in scaled space

## ---------- Back-transform each draw to original scale ----------
p_mu <- mean(p);   p_sd <- sd(p)
y_mu <- mean(y); y_sd <- sd(y)
pc_mu <- mean(pc); pc_sd <- sd(pc)
ps_mu <- mean(ps); ps_sd <- sd(ps)
alphaDem_draws <- cbind(
  alpha0  = postDem_s[,1] - postDem_s[,2] * (p_mu/p_sd) - postDem_s[,3] * (y_mu/y_sd)
  - postDem_s[,4] * (pc_mu/pc_sd) - postDem_s[,5] * (pc_mu/pc_sd),
  alphap  = postDem_s[,2] / p_sd,
  alphay = postDem_s[,3] / y_sd,
  alphapc = postDem_s[,4] / pc_sd,
  alphaps = postDem_s[,5] / ps_sd
)

## ---------- Summaries ----------
summ <- function(x) c(mean=mean(x), sd=sd(x), quantile(x, c(.025,.25,.5,.75,.975)))
cat("\nPosterior (original scale):\n")
print(rbind(
  alpha0  = summ(alphaDem_draws[,"alpha0"]),
  alphap  = summ(alphaDem_draws[,"alphap"]),
  alphay = summ(alphaDem_draws[,"alphay"]),
  alphapc = summ(alphaDem_draws[,"alphapc"]),
  alphaps = summ(alphaDem_draws[,"alphaps"])
))
cat("\nTrue (original scale):\n")
print(c(alpha0=B1, alphap=B2, alphay=B3, alphpc=B4, alphaps=B5))
ElastDemPrice <- alphaDem_draws[,2]
tax <- 0.1
CausalEffect <- (ElastSupPrice*ElastDemPrice)*log(1+tax)/(ElastSupPrice-ElastDemPrice) 
popCauEff <- (A2 * B2)/(A2 - B2) * log(1 + tax)

library(ggplot2)

post_draws <- as.numeric(CausalEffect)

dens <- density(post_draws, adjust = 1)
df   <- data.frame(x = dens$x, y = dens$y)
ci   <- quantile(post_draws, c(0.025, 0.975))
pm   <- mean(post_draws)

lines_df <- data.frame(
  x    = c(pm, as.numeric(popCauEff)),
  name = c("Posterior mean", "Population value")
)

ggplot(df, aes(x, y)) +
  geom_area(data = subset(df, x >= ci[1] & x <= ci[2]), alpha = 0.15) +
  geom_line(linewidth = 1) +
  geom_vline(data = lines_df,
             aes(xintercept = x, color = name, linetype = name),
             linewidth = 0.8) +
  labs(
    x = expression(tau), y = "Density",
    title = "Posterior of the causal effect",
    subtitle = sprintf("95%% CrI: [%.3f, %.3f]; mean = %.3f", ci[1], ci[2], pm),
    color = NULL, linetype = NULL
  ) +
  scale_linetype_manual(values = c("Posterior mean" = "solid",
                                   "Population value" = "dashed")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

############# Instrumental Variable Quantile Regression: Simulation ##############
rm(list = ls()); set.seed(10101)
library(quantreg)
### Simulation ####
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

######### Doubly Robust Bayesian: Simulation with misspecification #########
rm(list = ls()); set.seed(123)

## simulate covariates
n  <- 2000
X1 <- rnorm(n)                 # continuous
X2 <- rbinom(n, 1, 0.5)        # binary
X3 <- runif(n, -1, 1)          # continuous
Z <- cbind(X1, X2, X3)
# True propensity and outcome
logit <- function(t) 1/(1+exp(-t))

true_pi  <- function(X){ logit(-0.3 + 0.8*X[,1] -0.5*X[,2] + 0.7*X[,3] + 0.7*X[,1]*X[,3] + 0.5*sin(X[,1])) }
true_eta <- function(d,X){ -0.2 + 0.6*d + 0.5*X[,1] -0.4*X[,2] + 0.3*X[,3] + 0.6*(X[,1]^2) - 0.5*X[,1]*X[,2] }
true_p   <- function(d,X){ logit(true_eta(d,X)) }

ATE  <- mean(true_p(1,Z) - true_p(0,Z))

pi <- true_pi(Z)
D  <- rbinom(n,1,pi)
P  <- true_p(D,Z)
Y  <- rbinom(n,1,P)

dat <- data.frame(Y, D, X1, X2 = factor(X2), X3)

## fit parametric logit for treatment (propensity)
ps_fit  <- glm(D ~ X1 + X2 + X3 + X1:X3 + I(sin(X1)), data = dat, family = binomial(link = "logit"))
ps_hat  <- predict(ps_fit, type = "response")  # \hat{\pi}(X)

## fit parametric logit for binary outcome (with treatment)
out_fit <- glm(Y ~ D + X1 + X2 + X3 + I(X1^2) + X1:X2, data = dat, family = binomial(link = "logit"))
p_hat   <- predict(out_fit, type = "response") # \hat{p}(D, X)

summary(ps_fit)
summary(out_fit)

# potential-outcome predictions
mu1_hat <- predict(out_fit, newdata = transform(dat, D = 1), type = "response")
mu0_hat <- predict(out_fit, newdata = transform(dat, D = 0), type = "response")

ATE_hat <- mean(mu1_hat - mu0_hat)                 # risk difference (ATE)

# Standard error using bootstrap
set.seed(123)
B <- 500
n <- nrow(dat)

ATE_b <- replicate(B, {
  idx <- sample.int(n, n, replace = TRUE)
  db  <- dat[idx, ]
  # skip degenerate resamples with one treatment level
  if (length(unique(db$D)) < 2) return(NA_real_)
  fit <- glm(Y ~ D + X1 + X2 + X3 + I(X1^2) + X1:X2, data = db, family = binomial("logit"))
  p1  <- predict(fit, newdata = transform(db, D = 1), type = "response")
  p0  <- predict(fit, newdata = transform(db, D = 0), type = "response")
  mean(p1 - p0)
})

SE_boot <- sd(na.omit(ATE_b))
c(SE_boot = SE_boot)

######### Doubly Robust Bayesian: Implementation #########

library(glmnet)
# OK specification propensity score
# Consistent if the PS model is correctly specified (even if the outcome model is wrong). Requires positivity and often benefits from trimming/extreme-weight control.
Z <- cbind(X1, X2, X3, X1*X3, sin(X1))
p <- dim(Z)[2]
covars <- c("X1","X2","X3", "X1X3", "sinX1")

# # OK specification outcome regression
# # Consistent if the outcome model is correctly specified (even if the PS is wrong).
# Z <- cbind(X1, X2, X3, X1*X2, X1^2)
# p <- dim(Z)[2]
# covars <- c("X1","X2","X3", "X1X2", "X12")

Z.df <- as.data.frame(Z) # covariates as a dataframe
cvfit <- cv.glmnet(Z, D, family = "binomial", alpha=1, nfolds = 10) # lasso,min
ps.coef <- coef(cvfit, s = "lambda.min")[0:p+1]
ps_hat <- predict(cvfit, newx = Z, s = "lambda.min",type = "response")
# Avoid 0/1
ps_hat <- pmin(pmax(ps_hat, 1e-6), 1 - 1e-6)

# Riesz representer gamma(d,x) = d/\pi(x) - (1-d)/(1-\pi(x))
riesz_fun <- function(d, ps) d/ps - (1 - d)/(1 - ps)

# -------------------------------
# GP helpers
# -------------------------------
# Build a GP model with (optionally) an additive linear kernel on gamma
make_gp <- function(X_train, y_bin, use_correction = FALSE, gamma_scaled = NULL,
                    approx_method = gplite::approx_laplace(),
                    init_lscale = 0.3) {
  stopifnot(length(y_bin) == nrow(X_train))
  # Base SE kernel (isotropic over all columns provided)
  # Initializing the covariance function using squared exponential kernel
  cf_se <- gplite::cf_sexp(vars = colnames(X_train), lscale = init_lscale, magn = 1, normalize = TRUE)
  cfs <- list(cf_se)
  if (use_correction) {
    stopifnot(!is.null(gamma_scaled))
    # Add linear kernel on the last column (gamma_scaled)
    # Kernel: k_PS(i,j) = (gamma_scaled_i)*(gamma_scaled_j)
    cf_lin <- gplite::cf_lin(vars = "gamma", magn = 1, normalize = FALSE)
    cfs <- list(cf_se, cf_lin)
  }
  # Initializes a GP model with given covariance function(s) and likelihood
  gp <- gplite::gp_init(cfs = cfs, lik = gplite::lik_bernoulli(), approx = approx_method)
  # Optimizing hyperparameters
  gp <- gplite::gp_optim(gp, x = X_train, y = y_bin, maxiter = 1000, restarts = 3, tol = 1e-05)
  return(gp)
}

# Draw posterior of m(Z,d) at test inputs using GP classification
# Returns a draws x n_test matrix on probability scale (logistic link applied)
posterior_draws_prob <- function(gp, X_test, ndraws = 5000) {
  # NOTE: gplite::gp_draw returns an N_test x draws matrix.
  # We want draws x N_test to align with downstream code (weights W: N_post x n_eff).
  out <- gplite::gp_draw(gp, xnew = X_test, draws = ndraws,
                         transform = TRUE, target = FALSE, jitter = 1e-6)
  t(out)
}

# -------------------------------
# Main routine for a given trimming threshold t in {0.10, 0.05, 0.01}
# -------------------------------
run_trim <- function(t_trim = 0.10, N_post = 5000, h_r = 1/2, cor_size = 1) {
  # Trim by estimated PS
  keep <- which(ps_hat >= t_trim & ps_hat <= (1 - t_trim))
  n_eff <- length(keep)
  if (n_eff < 50) stop("Too few observations after trimming.")
  
  Yk  <- Y[keep]
  Dk  <- D[keep]
  Zk  <- Z[keep, , drop = FALSE]
  psk <- ps_hat[keep]
  
  # Riesz representer and sigma_n choice
  gamma_k <- riesz_fun(Dk, psk)
  # sigma_n ≍ log(n) / ( n^{1/2} * mean|gamma| )  (scale-matched)
  sigma_n <- cor_size * (log(n_eff) / ( (n_eff)^(h_r) * mean(abs(gamma_k)) ))
  gamma_scaled <- sigma_n * gamma_k
  
  # Training inputs for GP: [Z, D] and, if corrected, an extra column 'gamma'
  X_train_nc <- data.frame(Zk)
  X_train_nc$D <- Dk
  colnames(X_train_nc) <- c(covars, "D")
  
  X_train_c <- X_train_nc
  # Note that in the paper says that the adjusted prior has the factor \lambda*\hat{\gamma},
  # where \lambda \sim N(0,\sigma_n^2). There is multiplication by \lambda in the code
  # because this is integrated out. This approach is equivalent what we see in the paper
  X_train_c$gamma <- gamma_scaled
  
  # Test inputs: for each i, (Z_i, d=0) and (Z_i, d=1)
  X_t0 <- X_train_nc; X_t0$D <- 0
  X_t1 <- X_train_nc; X_t1$D <- 1
  X_t0c <- X_t0; X_t1c <- X_t1
  X_t0c$gamma <- sigma_n * riesz_fun(0, psk)   # = -sigma_n/(1-ps)
  X_t1c$gamma <- sigma_n * riesz_fun(1, psk)   # =  sigma_n/ps
  
  # ---------------- GP WITHOUT prior correction ----------------
  gp_nc <- make_gp(X_train = X_train_nc, y_bin = Yk, use_correction = FALSE)
  # joint draws for [m(Z,0); m(Z,1)]
  M_nc_0 <- posterior_draws_prob(gp_nc, X_t0, ndraws = N_post)  # N_post x n_eff
  M_nc_1 <- posterior_draws_prob(gp_nc, X_t1, ndraws = N_post)
  # observed arm probabilities
  M_nc_obs <- ifelse(matrix(Dk, nrow = N_post, ncol = n_eff, byrow = TRUE) == 1, M_nc_1, M_nc_0)
  
  # ---------------- GP WITH prior correction -------------------
  gp_c <- make_gp(X_train = X_train_c, y_bin = Yk, use_correction = TRUE, gamma_scaled = gamma_scaled)
  M_c_0 <- posterior_draws_prob(gp_c, X_t0c, ndraws = N_post)
  M_c_1 <- posterior_draws_prob(gp_c, X_t1c, ndraws = N_post)
  M_c_obs <- ifelse(matrix(Dk, nrow = N_post, ncol = n_eff, byrow = TRUE) == 1, M_c_1, M_c_0)
  
  # ---------------- Bayesian bootstrap weights -----------------
  W <- matrix(rexp(N_post * n_eff, rate = 1), nrow = N_post)
  W <- W / rowSums(W)
  
  # ---------------- ATEs ----------------
  # (1) Unadjusted Bayes
  Ate_nc_draws <- rowSums((M_nc_1 - M_nc_0) * W)
  # (2) Prior-adjusted Bayes
  Ate_c_draws  <- rowSums((M_c_1 - M_c_0) * W)
  # (3) DR Bayes (posterior recentering)
  # Pre-centering using UNcorrected GP posterior means
  mu_nc_diff <- colMeans(M_nc_1 - M_nc_0)                 # length n_eff
  mu_nc_obs  <- colMeans(M_nc_obs)                        # length n_eff
  ATE_dr_pre <- mean(mu_nc_diff + gamma_k * (Yk - mu_nc_obs))
  # Recenter draws using corrected GP draws
  DR_rec_1 <- (matrix(gamma_k, nrow = N_post, ncol = n_eff, byrow = TRUE)) * (matrix(Yk, nrow = N_post, ncol = n_eff, byrow = TRUE) - M_c_obs)
  # The first component is \tau_{\eta}^s in Algorithm 1
  # The second component is \hat{m} a scalar. This has positive sign because
  # minus x minus, the bias is with minus, and then, this component enters with minus in the bias term
  # The third component is \gamma_k x (y-m(d,x)), m(d,x) is with prior correction
  # The fourth component is m(1,x)-m(0,x) also with correction in the prior
  Ate_drb_draws <- rowSums((M_c_1 - M_c_0) * W) + ATE_dr_pre - rowSums(DR_rec_1) / n_eff - rowSums(M_c_1 - M_c_0) / n_eff
  
  qfun <- function(x) {
    qs <- quantile(x, c(0.025, 0.975))
    c(
      mean   = mean(x),
      median = median(x),
      lo     = qs[1],
      hi     = qs[2],
      len    = diff(qs)
    )
  }
  
  vals <- list(
    Bayes      = qfun(Ate_nc_draws),
    `PA Bayes` = qfun(Ate_c_draws),
    `DR Bayes` = qfun(Ate_drb_draws)
  )
  out_df <- as.data.frame(do.call(rbind, vals), check.names = FALSE)
  # # Guarantee expected columns exist
  # needed <- c("mean","lo","hi","len","median")
  # for (nm in needed) if (!nm %in% colnames(out_df)) out_df[[nm]] <- NA_real_
  # # Order the summary columns
  # out_df <- out_df[, needed, drop = FALSE]
  # # Add n_eff as a separate column (used only for the header line while printing)
  # out_df$n_eff <- n_eff
  list(t_trim = t_trim, results = out_df)
}

N_post <- 5000
h_r    <- 1/2
cor_size <- 1
trim <- 0.05
res_list <- run_trim(t_trim = trim, N_post = N_post, h_r = h_r, cor_size = cor_size)
print(res_list$results)
ATE; ATE_hat; sd(ATE_b)

