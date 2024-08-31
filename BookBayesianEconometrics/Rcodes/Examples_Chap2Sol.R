# Examples in chapter 2


########################## The math test ########################## 

N <- 50 # Sample size
y_bar <- 102 # Sample mean 
s2 <- 10 # Sample variance
alpha <- N - 1
serror <- (s2/N)^0.5 
y.H0 <- c(100, 100.5, 101, 101.5, 102)
test <- (y.H0 - y_bar)/serror
pval <- 2*pt(test, alpha)
pval
# p-values
PO01 <- (gamma(N/2)*((N-1)*serror^2)^(-0.5)*
           (1+test^2/alpha)^(-N/2))/(gamma(1/2)*gamma((N-1)/2))
PO01/(1+PO01)
# Posterior model probability of the null hypothesis.