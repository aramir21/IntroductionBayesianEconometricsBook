########################## Simulation exercise: AR(2) model ########################## 
rm(list = ls())
set.seed(010101)
T <- 150
mu <- 0.007; phi1 <- 0.5; phi2 <- 0.3; sig <- 0.035
Ey <- mu/(1-phi1-phi2); Sigy <- sig*((1-phi2)/(1-phi2-phi1^2-phi2*phi1^2-phi2^2+phi2^3))^0.5 
y <- rnorm(T, mean = Ey, sd = Sigy)
e <- rnorm(T, mean = 0, sd = sig)
for(t in 3:T){
  y[t] <- mu + phi1*y[t-1] + phi2*y[t-2] + e[t]
}
mean(y); sd(y)
y <- ts(y, start=c(1920, 1), frequency=1)
plot(y)
iter <- 1000; burnin <- 500
library(bayesforecast)
sf1 <- bayesforecast::stan_sarima(y,order = c(2, 0, 0), prior_mu0 = normal(0, 1),
                                  prior_ar = normal(0, 1), prior_sigma0 = inverse.gamma(0.01, 0.01),
                  seasonal = c(0, 0, 0), iter = iter, warmup = burnin, chains = 1)
Postmu <- sf1[["stanfit"]]@sim[["samples"]][[1]][["mu0"]][-c(1:burnin)]
Postsig <- sf1[["stanfit"]]@sim[["samples"]][[1]][["sigma0"]][-c(1:burnin)]
Postphi1 <- sf1[["stanfit"]]@sim[["samples"]][[1]][["ar0[1]"]][-c(1:burnin)]
Postphi2 <- sf1[["stanfit"]]@sim[["samples"]][[1]][["ar0[2]"]][-c(1:burnin)]
Postdraws <- cbind(Postmu, Postsig, Postphi1, Postphi2)
summary(coda::mcmc(Postdraws))
# print(sf1)
