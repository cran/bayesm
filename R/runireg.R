runireg=
function(y,X,betabar,A,nu,ssq){
#   
#    revision history:
#      p. rossi 1/11/05 -- fix error in sum of squares
#  
# purpose: 
#          draw from posterior for a univariate regression model
#          with natural conjugate prior
#
# arguments: 
#      y is n x 1
#      X is n x k
#      betabar is prior mean 
#      A is prior precision
#      nu, ssq are parameters of prior on sigmasq
# output:
#      list of beta, sigmasq draws
#         beta is k x 1 vector of coefficients
# model:
#      Y=Xbeta+e  var(e_i) = sigmasq
#      priors:  beta| sigmasq ~ N(betabar,sigmasq*A^-1)
#                sigmasq ~ (nu*ssq)/chisq_nu
n=length(y)
k=ncol(X)
#
# first draw Sigma
#
RA=chol(A)
W=rbind(X,RA)
z=c(y,as.vector(RA%*%betabar))
IR=backsolve(chol(crossprod(W)),diag(k))
#      W'W=R'R ;  (W'W)^-1 = IR IR'  -- this is UL decomp
btilde=crossprod(t(IR))%*%crossprod(W,z)
res=z-W%*%btilde
s=t(res)%*%res
#
# first draw Sigma
#
#
sigmasq=(nu*ssq + s)/rchisq(1,nu+n)
#
# now draw beta given Sigma
#	
beta = btilde + as.vector(sqrt(sigmasq))*IR%*%rnorm(k)
list(beta=beta,sigmasq=sigmasq)
}
