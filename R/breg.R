breg=
function(y,X,betabar,A) 
{
#
# P.Rossi 12/04
#  revision history:
#    P. Rossi 3/27/05 -- changed to augment strategy
#
# Purpose: draw from posterior for linear regression, sigmasq=1.0
#
# Output:  draw from posterior
# 
# Model: y = Xbeta + e  e ~ N(0,I)
#
# Prior:  beta ~ N(betabar,A^-1)
#
k=length(betabar)
RA=chol(A)
W=rbind(X,RA)
z=c(y,as.vector(RA%*%betabar))
IR=backsolve(chol(crossprod(W)),diag(k))
#      W'W=R'R ;  (W'W)^-1 = IR IR'  -- this is UL decomp
return(crossprod(t(IR))%*%crossprod(W,z)+IR%*%rnorm(k))

}
