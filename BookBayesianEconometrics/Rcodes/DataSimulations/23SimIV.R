set.seed(66)
simIV = function(delta,beta,Sigma,n,z,w,gamma) {
  eps = matrix(rnorm(2*n),ncol=2) %*% chol(Sigma)
  x = z %*% delta + eps[,1]; y = beta*x + eps[,2] + w%*%gamma
  list(x=as.vector(x),y=as.vector(y)) }
n = 1000 ; p=1 # number of instruments
z = cbind(rep(1,n),matrix(runif(n*p),ncol=p))
w = matrix(1,n,1)
rho=.8
Sigma = matrix(c(1,rho,rho,1),ncol=2)
delta = c(1,4); beta = .5; gamma = c(1)
simiv = simIV(delta,beta,Sigma,n,z,w,gamma)
Dat<-cbind(simiv$y,simiv$x,z[,2])
colnames(Dat)<-c("y","x","z")
write.csv(Dat,file="23SimIV.csv")


