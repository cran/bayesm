rhierLinearModel=
function(Data,Prior,Mcmc)
{
#
# Revision History
#     1/17/05  P. Rossi
#
# Purpose:
#   run hiearchical regression model
#
# Arguments:
#   Data list of regdata,Z 
#     regdata is a list of lists each list with members y, X
#        e.g. regdata[[i]]=list(y=y,X=X)
#     X has nvar columns
#     Z is nreg=length(regdata) x nz
#   Prior list of prior hyperparameters
#     Deltabar,A, nu.e,ssq,nu,V
#          note: ssq is a nreg x 1 vector!
#   Mcmc
#     list of Mcmc parameters
#     R is number of draws
#     keep is thining parameter -- keep every keepth draw
#
# Output: 
#   list of 
#   betadraw -- nreg x nvar x R/keep array of individual regression betas
#   taudraw -- R/keep x nreg  array of error variances for each regression
#   Deltadraw -- R/keep x nz x nvar array of Delta draws
#   Vbetadraw -- R/keep x nvar*nvar array of Vbeta draws
#
# Model:
# nreg regression equations 
#        y_i = X_ibeta_i + epsilon_i  
#        epsilon_i ~ N(0,tau_i)
#             nvar X vars in each equation
#
# Priors:
#        tau_i ~ nu.e*ssq_i/chisq(nu.e)  tau_i is the variance of epsilon_i
#        beta_i ~ N(ZDelta[i,],V_beta)
#               Note:  ZDelta is the matrix Z * Delta; [i,] refers to ith row of this product!
#
#          vec(Delta) | V_beta ~ N(vec(Deltabar),Vbeta (x) A^-1)
#          V_beta ~ IW(nu,V)  or V_beta^-1 ~ W(nu,V^-1)
#              Delta, Deltabar are nz x nvar
#              A is nz x nz
#              Vbeta is nvar x nvar
#        
#          NOTE: if you don't have any z vars, set Z=iota (nreg x 1) 
#
#
#  create needed functions
#
#------------------------------------------------------------------------------
append=function(l) { l=c(l,list(XpX=crossprod(l$X),Xpy=crossprod(l$X,l$y)))}
#
getvar=function(l) { var(l$y)}
#
runiregG=
function(y,X,XpX,Xpy,sigmasq,A,betabar,nu,ssq){
# 
# Purpose:
#   perform one Gibbs iteration for Univ Regression Model
#   only does one iteration so can be used in rhierLinearModel
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
n=length(y)
k=ncol(XpX)
sigmasq=as.vector(sigmasq)
#
#     first draw beta | sigmasq
#
  IR=backsolve(chol(XpX/sigmasq+A),diag(k))
  btilde=crossprod(t(IR))%*%(Xpy/sigmasq+A%*%betabar)
  beta = btilde + IR%*%rnorm(k)
#
#    now draw sigmasq | beta
#
  res=y-X%*%beta
  s=t(res)%*%res
  sigmasq=(nu*ssq + s)/rchisq(1,nu+n)

list(betadraw=beta,sigmasqdraw=sigmasq)
}

#------------------------------------------------------------------------------
#

#
# check arguments
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata and Z")}
    if(is.null(Data$regdata)) {pandterm("Requires Data element regdata")}
    regdata=Data$regdata
    if(is.null(Data$Z)) {pandterm("Requires Data element Z")}
    Z=Data$Z
nz=ncol(Z)
nreg=length(regdata)
#
# check data for validity
#
dimfun=function(l) {c(length(l$y),dim(l$X))}
dims=sapply(regdata,dimfun)
dims=t(dims)
nvar=quantile(dims[,3],prob=.5)

for (i in 1:nreg) 
{
   if(dims[i,1] != dims[i,2]  || dims[i,3] !=nvar) 
      {pandterm(paste("Bad Data dimensions for unit ",i," dims(y,X) =",dims[i,]))}
}
#
# check for Prior
#
if(missing(Prior))
   { Deltabar=matrix(rep(0,nz*nvar),ncol=nvar); A=diag(c(rep(.01,nz)));
     nu.e=3; ssq=sapply(regdata,getvar) ; nu=nvar+3 ; V= nu*diag(nvar)}
else
   {
    if(is.null(Prior$Deltabar)) {Deltabar=matrix(rep(0,nz*nvar),ncol=nvar)} 
       else {Deltabar=Prior$Deltabar}
    if(is.null(Prior$A)) {A=diag(c(rep(.01,nz)))} 
       else {A=Prior$A}
    if(is.null(Prior$nu.e)) {nu.e=3} 
       else {nu.e=Prior$nu.e}
    if(is.null(Prior$ssq)) {ssq=sapply(regdata,getvar)} 
       else {ssq=Prior$ssq}
    if(is.null(Prior$nu)) {nu=nvar+3} 
       else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(nvar)} 
       else {V=Prior$V}
   }
#
# check dimensions of Priors
#
if(ncol(A) != nrow(A) || ncol(A) != nz || nrow(A) != nz) 
   {pandterm(paste("bad dimensions for A",dim(A)))}
if(nrow(Deltabar) != nz || ncol(Deltabar) != nvar)
   {pandterm(paste("bad dimensions for Deltabar ",dim(Deltabar)))}
if(length(ssq) != nreg) {pandterm(paste("bad length for ssq ",length(ssq)))}
if(ncol(V) != nvar || nrow(V) != nvar) {pandterm(paste("bad dimensions for V ",dim(V)))}
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
cat("Starting Gibbs Sampler for Linear Hierarchical Model",fill=TRUE)
cat("   ",nreg," Regressions",fill=TRUE)
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("Deltabar",fill=TRUE)
print(Deltabar)
cat("A",fill=TRUE)
print(A)
cat("nu.e (d.f. parm for regression error variances)= ",nu.e,fill=TRUE)
cat("Vbeta ~ IW(nu,V)",fill=TRUE)
cat("nu = ",nu,fill=TRUE)
cat("V ",fill=TRUE)
print(V)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep,fill=TRUE)
cat(" ",fill=TRUE)
#
#  allocate space for the draws and set initial values of Vbeta and Delta
#
Vbetadraw=matrix(double(floor(R/keep)*nvar*nvar),ncol=nvar*nvar)
Deltadraw=matrix(double(floor(R/keep)*nz*nvar),ncol=nz*nvar)
taudraw=matrix(double(floor(R/keep)*nreg),ncol=nreg)
betadraw=array(double(floor(R/keep)*nreg*nvar),dim=c(nreg,nvar,floor(R/keep)))

tau=double(nreg)
Delta=c(rep(0,nz*nvar))
Vbeta=as.vector(diag(nvar))
betas=matrix(double(nreg*nvar),ncol=nvar)
#
#  set up fixed parms for the draw of Vbeta,Delta
#
#  note: in the notation of the MVR  Y =    X      B  
#                                  n x m  n x k  k x m
#                           "n" = nreg
#                           "m" = nvar
#                           "k" = nz
#			general model: Beta = Z Delta + U
#
Fparm=init.rmultiregfp(Z,A,Deltabar,nu,V)
#
#       Create XpX elements of regdata and initialize tau
#
regdata=lapply(regdata,append)

tau=sapply(regdata,getvar)
#
#	start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
fsh()

for(rep in 1:R)
{
   Abeta=chol2inv(chol(matrix(Vbeta,ncol=nvar)))
   betabar=Z%*%matrix(Delta,ncol=nvar)
#
#       loop over all regressions
#
   for (reg in 1:nreg) 
   {
      regout=runiregG(regdata[[reg]]$y,regdata[[reg]]$X,regdata[[reg]]$XpX,
                regdata[[reg]]$Xpy,tau[reg],Abeta,betabar[reg,],nu.e,ssq[reg])
      betas[reg,]=regout$betadraw
      tau[reg]=regout$sigmasqdraw
   }
#
#          draw Vbeta, Delta | {beta_i}
#
   rmregout=rmultiregfp(betas,Z,Fparm)
   Vbeta=as.vector(rmregout$Sigma)
   Delta=as.vector(rmregout$B)
#
#       print time to completion and draw # every 100th draw
#
  if(rep%%100 == 0)
    {ctime=proc.time()[3]
    timetoend=((ctime-itime)/rep)*(R-rep)
    cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
    fsh()}


  if(rep%%keep == 0) 
    {mkeep=rep/keep
     Vbetadraw[mkeep,]=Vbeta
     Deltadraw[mkeep,]=Delta
     taudraw[mkeep,]=tau
     betadraw[,,mkeep]=betas}

}
ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')

return(list(Vbetadraw=Vbetadraw,Deltadraw=Deltadraw,betadraw=betadraw,taudraw=taudraw))
}
