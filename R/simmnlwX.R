simmnlwX=
function(n,X,beta) 
{
#
#     p. rossi 2004
#
# Purpose:  simulate from MNL model conditional on X matrix
#
# Arguments:
#   j is number of alternatives
#   X is full matrix, including intercepts
#   beta is true parm value
#
# Output:
#   list of X  (note: we include full set of intercepts and 2 unif(-1,1) X vars)
#   y  (indicator of choice-- 1, ...,j
#   prob is a n x j matrix of choice probs
#
#   note: first choice alternative has intercept set to zero
#
#
k=length(beta)
# now construct probabilities
Xbeta=X%*%beta
j=nrow(Xbeta)/n
Xbeta=matrix(Xbeta,byrow=T,ncol=j)
Prob=exp(Xbeta)
iota=c(rep(1,j))
denom=Prob%*%iota
Prob=Prob/as.vector(denom)
# draw y
y=vector("double",n)
ind=1:j
for (i in 1:n) 
        {
        yvec=rmultinom(1,1,Prob[i,])
        y[i]=ind%*%yvec
        }

list(y=y,X=X,beta=beta,prob=Prob)
}
