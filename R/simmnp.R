simmnp=
function(X,p,n,beta,sigma) {
#
# simulate from MNP:  w_i=X_i*beta + e   e ~ MVN(0,Sigma)
#                     y_i = j  if max(w_i) = w_ij
#			       if max(w_i) < 0  y_i = p
#	Sigma is p-1 x p-1
#       X is n(p-1) x k
#       beta is k x 1
#
#   create functions needed
#
indmax=function(x) {
ind=1:length(x)
ind[x==max(x)]
}

Xbeta=X%*%beta
w=as.vector(crossprod(chol(sigma),matrix(rnorm((p-1)*n),ncol=n)))+ Xbeta
w=matrix(w,ncol=(p-1),byrow=T)
maxw=apply(w,1,max)
y=apply(w,1,indmax)
y=ifelse(maxw < 0,p,y)

return(list(y=y,X=X,beta=beta,sigma=sigma))
}

