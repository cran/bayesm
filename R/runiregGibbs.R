runiregGibbs=
function(Data,Prior,Mcmc)
{
# 
# revision history:
#          P. Rossi 1/17/05
# Purpose:
#   perform Gibbs iterations for Univ Regression Model using
#     prior with beta, sigma-sq indep
# 
# Arguments:
#   Data -- list of data 
#           y,X
#   Prior -- list of prior hyperparameters
#     betabar,A      prior mean, prior precision
#     nu, ssq        prior on sigmasq
#   Mcmc -- list of MCMC parms
#     sigmasq=initial value for sigmasq
#     R number of draws
#     keep -- thinning parameter
# 
# Output: 
#   list of beta, sigmasq
#
# Model:
#   y = Xbeta + e  e ~N(0,sigmasq)
#          y is n x 1
#          X is n x k
#          beta is k x 1 vector of coefficients
#
# Priors:  beta ~ N(betabar,A^-1)
#          sigmasq ~ (nu*ssq)/chisq_nu
# 
#
# check arguments
#
pandterm=function(message) {stop(message,call.=FALSE)}
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
if(nobs != nrow(X) ) {pandterm("length(y) ne nrow(X)")}
#
# check for Prior
#
if(missing(Prior))
   { betabar=c(rep(0,nvar)); A=.01*diag(nvar)}
else
   {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))} 
       else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=.01*diag(nvar)} 
       else {A=Prior$A}
    if(is.null(Prior$nu)) {nu=3}
       else {nu=Prior$nu}
    if(is.null(Prior$ssq)) {ssq=var(y)}
       else {ssq=Prior$ssq}
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
    if(is.null(Mcmc$sigmasq)) {sigmasq=var(y)} else {sigmasq=Mcmc$sigmasq}
   }
#
# print out problem
#
cat(" ", fill=TRUE)
cat("Starting Gibbs Sampler for Univariate Regression Model",fill=TRUE)
cat("  with ",nobs," observations",fill=TRUE)
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat("nu = ",nu," ssq= ",ssq,fill=TRUE)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep,fill=TRUE)
cat(" ",fill=TRUE)

sigmasqdraw=double(floor(Mcmc$R/keep))
betadraw=matrix(double(floor(Mcmc$R*nvar/keep)),ncol=nvar)
XpX=crossprod(X)
Xpy=crossprod(X,y)
sigmasq=as.vector(sigmasq)

itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
fsh()

for (rep in 1:Mcmc$R)
{
#
#     first draw beta | sigmasq
#
  IR=backsolve(chol(XpX/sigmasq+A),diag(nvar))
  btilde=crossprod(t(IR))%*%(Xpy/sigmasq+A%*%betabar)
  beta = btilde + IR%*%rnorm(nvar)
#
#    now draw sigmasq | beta
#
  res=y-X%*%beta
  s=t(res)%*%res
  sigmasq=(nu*ssq + s)/rchisq(1,nu+nobs)
  sigmasq=as.vector(sigmasq)
#
#       print time to completion and draw # every 100th draw
#
  if(rep%%100 == 0)
    {ctime=proc.time()[3]
    timetoend=((ctime-itime)/rep)*(R-rep)
    cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
    fsh()}

  if(rep%%keep == 0) 
    {mkeep=rep/keep; betadraw[mkeep,]=beta; sigmasqdraw[mkeep]=sigmasq}
}
ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')


return(list(betadraw=betadraw,sigmasqdraw=sigmasqdraw))
}
