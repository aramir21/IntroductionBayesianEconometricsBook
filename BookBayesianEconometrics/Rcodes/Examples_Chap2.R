# Examples in chapter 2


########################## The likelihood principle ########################## 

success <- 9 # Number of observed success in n trials
n <- 12 # Number of trials
siglevel <- sum(sapply(9:n, function(y) dbinom(y, n, 0.5)))
paste("Significance level from binomial model is", sep = " ", round(siglevel, 4))

success <- 3 # Number of target success (tails)
failures <- 9 # Number of failures
siglevel <- 1 - pnbinom((failures - 1), success, 0.5)
paste("Significance level from negative binomial model is", sep = " ", round(siglevel, 4))

########################## A simple working example ##########################

N <- 50
# Sample size
y_bar <- 102 
#sample mean
s2 <- 10
#sample variance
alpha <- N - 1
serror <- (s2/N)^0.5 
LimInf <- y_bar - abs(qt(0.025, alpha)) * serror
# Lower bound
LimSup <- y_bar + abs(qt(0.025, alpha)) * serror
# Upper bound
paste("The 95% credible and confidence intervals are the same", "(", LimInf, LimSup, ").", "However, their interpretations are totally different", sep = " ")

y.cut <- 103
P <- 1-metRology::pt.scaled(y.cut, df = alpha, mean = y_bar, sd = serror)
paste("The probability that mu is greater than 103 is", P, sep = " ")


