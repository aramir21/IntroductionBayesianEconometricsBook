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


