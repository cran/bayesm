lndIChisq=
function(nu,ssq,x)
{
#
# P. Rossi 12/04
#
# Purpose: evaluate log-density of scaled Inverse Chi-sq
#  density of r.var. Z=nu*ssq/chisq(nu)
#
return(-lgamma(nu/2)+(nu/2)*log((nu*ssq)/2)-((nu/2)+1)*log(x)-(nu*ssq)/(2*x))
}
