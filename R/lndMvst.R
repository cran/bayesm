lndMvst=
function(x,nu,mu,rooti,NORMC=FALSE)
{
#
# modified by Rossi 12/2005 to include normalizing constant
#
# function to evaluate log of MVstudent t density with nu df, mean mu,
# and with sigmai=rooti%*%t(rooti)   note: this is the UL decomp of sigmai not LU!
# rooti is in the inverse of upper triangular chol root of sigma
# or Sigma=root'root   root=inv(rooti)
#
dim=length(x)
if(NORMC) 
   {constant=(nu/2)*log(nu)+lgamma((nu+dim)/2)-(dim/2)*log(pi)-lgamma(nu/2)}
  else
   {constant=0}
z=as.vector(t(rooti)%*%(x-mu))
return(constant -((dim+nu)/2)*log(nu+z%*%z)+sum(log(diag(rooti))))
}
