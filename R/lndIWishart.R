lndIWishart=
function(nu,V,IW)
{
# 
# P. Rossi 12/04
#
# purpose: evaluate log-density of inverted Wishart
#    includes normalizing constant
#
# arguments:
#   nu is d. f. parm
#   V is location matrix
#   IW is the value at which the density should be evaluated
#
# output:
#   value of log density
#
# note: in this parameterization, E[IW]=V/(nu-k-1)
#
k=ncol(V)
Uiw=chol(IW)
lndetVd2=sum(log(diag(chol(V))))
lndetIWd2=sum(log(diag(Uiw)))
#
# first evaluate constant
#
const=((nu*k)/2)*log(2)+((k*(k-1))/4)*log(pi)
arg=(nu+1-c(1:k))/2
const=const+sum(lgamma(arg))
return(-const+nu*lndetVd2-(nu+k+1)*lndetIWd2-.5*sum(diag(V%*%chol2inv(Uiw))))
}
