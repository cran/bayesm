lndMvn=
function(x,mu,rooti)
{
#
# changed 12/05 by Rossi to include normalizing constant
#
# function to evaluate log of MV NOrmal density with  mean mu, var Sigma
# Sigma=t(root)%*%root   (root is upper tri cholesky root)
# Sigma^-1=rooti%*%t(rooti)   
# rooti is in the inverse of upper triangular chol root of sigma
#          note: this is the UL decomp of sigmai not LU!
#                Sigma=root'root   root=inv(rooti)
#
z=as.vector(t(rooti)%*%(x-mu))
return(  -(length(x)/2)*log(2*pi) -.5*(z%*%z) + sum(log(diag(rooti))))
}
