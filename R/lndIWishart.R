lndIWishart=
function(nu,S,IW)
{
# 
# P. Rossi 12/04
#
# purpose: evaluate log-density of inverted Wishart
#    includes normalizing constant
#
# arguments:
#   nu is d. f. parm
#   S is location matrix
#   IW is the value at which the density should be evaluated
#
# output:
#   value of log density
#
# note: in this parameterization, E[IW]=S/(nu-k-1)
#
k=ncol(S)
Uiw=chol(IW)
lndetSd2=sum(log(diag(chol(S))))
lndetIWd2=sum(log(diag(Uiw)))
#
# first evaluate constant
#
const=((nu*k)/2)*log(2)+((k*(k-1))/4)*log(pi)
arg=(nu+1-c(1:k))/2
const=const+sum(lgamma(arg))
return(-const+nu*lndetSd2-(nu+k+1)*lndetIWd2-.5*sum(diag(S%*%chol2inv(Uiw))))
}
