lndMvst=
function(x,nu,mu,rooti)
{
#
# function to evaluate log of MVstudent t density with nu df, mean mu,
# and with sigmai=rooti%*%t(rooti)   note: this is the UL decomp of sigmai not LU!
# rooti is in the inverse of upper triangular chol root of sigma
# or Sigma=root'root   root=inv(rooti)
#
z=as.vector(t(rooti)%*%(x-mu))
return(-((length(x)+nu)/2)*log(nu+z%*%z)+sum(log(diag(rooti))))
}
