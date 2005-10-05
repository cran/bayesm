llmnl= 
function(beta,y,X) 
{
#    p. rossi 2004
#    changed order of arguments to put beta first 9/05
#
# Purpose:evaluate log-like for MNL
#
# Arguments:
#   y is n vector with element = 1,...,j indicating which alt chosen
#   X is nj x k matrix of xvalues for each of j alt on each of n occasions
#   beta is k vector of coefs
# 
# Output: value of loglike
#
n=length(y)
j=nrow(X)/n
Xbeta=X%*%beta
Xbeta=matrix(Xbeta,byrow=T,ncol=j)
ind=cbind(c(1:n),y)
xby=Xbeta[ind]
Xbeta=exp(Xbeta)
iota=c(rep(1,j))
denom=log(Xbeta%*%iota)
return(sum(xby-denom))
}
