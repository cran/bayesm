lndMvn=
function(x,mu,rooti)
{
#
# function to evaluate log of MV NOrmal density with  mean mu, var Sigma
# and with sigma^-1=rooti%*%t(rooti)   
# rooti is in the inverse of upper triangular chol root of sigma
#          note: this is the UL decomp of sigmai not LU!
#                Sigma=root'root   root=inv(rooti)
#
z=as.vector(t(rooti)%*%(x-mu))
return(-.5*(z%*%z) + sum(log(diag(rooti))))
}
