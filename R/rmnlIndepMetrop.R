rmnlIndepMetrop=
function(Data,Prior,Mcmc)
{
#
# revision history:
#   p. rossi 1/05
#   2/9/05 fixed error in Metrop eval
#
# purpose: 
#   draw from posterior for MNL using Independence Metropolis
#
# Arguments:
#   Data - list of p,y,X  
#     p is number of alternatives
#     X is nobs*p x nvar matrix
#     y is nobs vector of values from 1 to p
#   Prior - list of A, betabar
#     A is nvar x nvar prior preci matrix
#     betabar is nvar x 1 prior mean
#   Mcmc
#     R is number of draws
#     keep is thinning parameter
#     nu degrees of freedom parameter for independence 
#        sampling density
#
# Output:
#   list of betadraws
#
# Model:   Pr(y=j) = exp(x_j'beta)/sum(exp(x_k'beta)
#
# Prior:   beta ~ N(betabar,A^-1)
#
# check arguments
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of p, y, X")}
    if(is.null(Data$X)) {pandterm("Requires Data element X")}
    X=Data$X
    if(is.null(Data$y)) {pandterm("Requires Data element y")}
    y=Data$y
    if(is.null(Data$p)) {pandterm("Requires Data element p")}
    m=Data$m
nvar=ncol(X)
nobs=length(y)
#
# check data for validity
#
if(length(y) != (nrow(X)/p) ) {pandterm("length(y) ne nrow(X)/p")}
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
    if(is.null(Mcmc$nu)) {nu=6} else {nu=Mcmc$nu}
   }
#
# print out problem
#
cat(" ", fill=TRUE)
cat("Starting Independence Metropolis Sampler for Multinomial Logit Model",fill=TRUE)
cat("  with ",p," alternatives",fill=TRUE)
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep," nu (df for st candidates) = ",nu,fill=TRUE)
cat(" ",fill=TRUE)

betadraw=matrix(double(floor(R/keep)*nvar),ncol=nvar)
#
# compute required quantities for indep candidates
#
beta=c(rep(0,nvar))
mle=optim(beta,llmnl,X=X,y=y,method="BFGS",hessian=TRUE,control=list(fnscale=-1))
beta=mle$par
betastar=mle$par
mhess=mnlHess(y,X,beta)
candcov=chol2inv(chol(mhess))
root=chol(candcov)
rooti=backsolve(root,diag(nvar))
priorcov=chol2inv(chol(A))
rootp=chol(priorcov)
rootpi=backsolve(rootp,diag(nvar))

#
#	start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
fsh()

oldlpost=llmnl(y,X,beta)+lndMvn(beta,betabar,rootpi)
oldlimp=lndMvst(beta,nu,betastar,rooti)
#       note: we don't need the determinants as they cancel in
#       computation of acceptance prob
naccept=0

for (rep in 1:R) 
{
   betac=rmvst(nu,betastar,root)
   clpost=llmnl(y,X,betac)+lndMvn(betac,betabar,rootpi)
   climp=lndMvst(betac,nu,betastar,rooti)
   ldiff=clpost+oldlimp-oldlpost-climp
   alpha=min(1,exp(ldiff))
   if(alpha < 1) {unif=runif(1)} else {unif=0}
   if (unif <= alpha)
      { beta=betac
        oldlpost=clpost
        oldlimp=climp
        naccept=naccept+1}
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
return(list(betadraw=betadraw,acceptr=naccept/R))
}
