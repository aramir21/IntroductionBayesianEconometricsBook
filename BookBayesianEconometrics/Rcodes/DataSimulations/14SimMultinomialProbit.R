simmnp = function(X, p, n, beta, sigma) {
  indmax = function(x) {which(max(x)==x)}
  Xbeta = X%*%beta
  w = as.vector(crossprod(chol(sigma),matrix(rnorm((p-1)*n),ncol=n))) + Xbeta
  w = matrix(w, ncol=(p-1), byrow=TRUE)
  maxw = apply(w, 1, max)
  y = apply(w, 1, indmax)
  y = ifelse(maxw < 0, p, y)
  return(list(y=y, X=X, beta=beta, sigma=sigma))
}
p = 3
n = 500
beta = c(-1,1,1,2)
Sigma = matrix(c(1, 0.5, 0.5, 1), ncol=2)
k = length(beta)
X1 = matrix(runif(n*p,min=0,max=2),ncol=p)
X2 = matrix(runif(n*p,min=0,max=2),ncol=p)
X = bayesm::createX(p,na=2,nd=NULL,Xa=cbind(X1,X2),Xd=NULL,DIFF=TRUE,base=p)
simout = simmnp(X,p,500,beta,Sigma)
SimMultProbNew<-cbind(simout$y,X1,X2)
colnames(SimMultProbNew)<-c("y","x11","x12","x13","x21","x22","x23")
write.csv(Dat,file="14SimMultProbmodelNew.csv", row.names=FALSE)

