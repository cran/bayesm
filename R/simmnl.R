simmnl=
function(p,n,beta) 
{
#
#   p. rossi 2004
#
# Purpose: simulate from MNL (including X values)
#
# Arguments:
#   p is number of alternatives
#   n is number of obs
#   beta is true parm value
#
# Output:
#   list of X  (note: we include full set of intercepts and 2 unif(-1,1) X vars)
#   y  (indicator of choice-- 1, ...,p
#   prob is a n x p matrix of choice probs
#
#   note: first choice alternative has intercept set to zero
#
k=length(beta)
x1=runif(n*p,min=-1,max=1)
x2=runif(n*p,min=-1,max=1)
I2=diag(rep(1,p-1))
zero=rep(0,p-1)
xadd=rbind(zero,I2)
for(i in 2:n) {
        xadd=rbind(xadd,zero,I2)
}

X=cbind(xadd,x1,x2)

# now construct probabilities
Xbeta=X%*%beta
p=nrow(Xbeta)/n
Xbeta=matrix(Xbeta,byrow=T,ncol=p)
Prob=exp(Xbeta)
iota=c(rep(1,p))
denom=Prob%*%iota
Prob=Prob/as.vector(denom)
# draw y
y=vector("double",n)
ind=1:p
for (i in 1:n) 
        {
        yvec=rmultinom(1,1,Prob[i,])
        y[i]=ind%*%yvec
        }

return(list(y=y,X=X,beta=beta,prob=Prob))
}
