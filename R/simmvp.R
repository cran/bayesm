simmvp=
function(X,p,n,beta,sigma)
{
#
# simulate from MVP:  w_i=X_i*beta + e_i   e ~ MVN(0,Sigma)
#                     y_ij = 1  if w_ij > 0
#
#	Sigma is p x p
#       X is np x k
#       beta is k x 1

#
Xbeta=X%*%beta
w=as.vector(crossprod(chol(sigma),matrix(rnorm(p*n),ncol=n)))+ Xbeta
y=ifelse(w<0,0,1)

return(list(y=y,X=X,beta=beta,sigma=sigma))
}
