rbprobitGibbs=
function(Data,Prior,Mcmc)
{
#
# revision history:
#   p. rossi 1/05
#
# purpose: 
#   draw from posterior for binary probit using Gibbs Sampler
#
# Arguments:
#   Data - list of X,y  
#     X is nobs x nvar, y is nobs vector of 0,1
#   Prior - list of A, betabar
#     A is nvar x nvar prior preci matrix
#     betabar is nvar x 1 prior mean
#   Mcmc
#     R is number of draws
#     keep is thinning parameter
#
# Output:
#   list of betadraws
#
# Model:   y = 1 if  w=Xbeta + e   > 0  e ~N(0,1)
#
# Prior:   beta ~ N(betabar,A^-1)
#
#
# ----------------------------------------------------------------------
# define functions needed
#
breg1=
function(root,X,y,Abetabar) 
{
#
#     p.rossi 12/04
#
# Purpose: draw from posterior for linear regression, sigmasq=1.0
# 
# Arguments:
#  root is chol((X'X+A)^-1)
#  Abetabar = A*betabar
#
# Output:  draw from posterior
# 
# Model: y = Xbeta + e  e ~ N(0,I)
#
# Prior:  beta ~ N(betabar,A^-1)
#
cov=crossprod(root,root)
betatilde=cov%*%(crossprod(X,y)+Abetabar)
betatilde+t(root)%*%rnorm(length(betatilde))
}

pandterm=function(message) {stop(message,call.=FALSE)}
#
# ----------------------------------------------------------------------
#
# check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of y and X")}
    if(is.null(Data$X)) {pandterm("Requires Data element X")}
    X=Data$X
    if(is.null(Data$y)) {pandterm("Requires Data element y")}
    y=Data$y
nvar=ncol(X)
nobs=length(y)
#
# check data for validity
#
if(length(y) != nrow(X) ) {pandterm("y and X not of same row dim")}
#
# check for Prior
#
if(missing(Prior))
   { betabar=c(rep(0,nvar)); A=diag(rep(.01,nvar))}
else
   {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))} 
       else {Deltabar=Prior$Deltabar}
    if(is.null(Prior$A)) {A=diag(rep(.01,nvar))} 
       else {A=Prior$A}
   }
#
# check dimensions of Priors
#
if(ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) 
   {pandterm(paste("bad dimensions for A",dim(A)))}
if(length(betabar) != nvar)
   {pandterm(paste("betabar wrong length, length= ",length(betabar)))}
#
# check MCMC argument
#
if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
else
   {
    if(is.null(Mcmc$R)) 
       {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
    }
#
# print out problem
#
cat(" ", fill=TRUE)
cat("Starting Gibbs Sampler for Binary Probit Model",fill=TRUE)
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep,fill=TRUE)
cat(" ",fill=TRUE)

betadraw=matrix(double(floor(R/keep)*nvar),ncol=nvar)
beta=c(rep(0,nvar))
sigma=c(rep(1,nrow(X)))
root=chol(chol2inv(chol((crossprod(X,X)+A))))
Abetabar=crossprod(A,betabar)
        a=ifelse(y == 0,-100, 0)
        b=ifelse(y == 0, 0, 100)
#
#	start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
fsh()

for (rep in 1:R) 
{
  # draw z given beta(i-1)
  mu=X%*%beta
  z=rtrun(mu,sigma,a,b)
  beta=breg1(root,X,z,Abetabar)
#
#       print time to completion and draw # every 100th draw
#
  if(rep%%100 == 0)
    {ctime=proc.time()[3]
    timetoend=((ctime-itime)/rep)*(R-rep)
    cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
    fsh()}

  if(rep%%keep == 0) 
    {mkeep=rep/keep; betadraw[mkeep,]=beta}
}
ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')
return(list(betadraw=betadraw))
}
