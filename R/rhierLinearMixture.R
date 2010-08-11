rhierLinearMixture=
function(Data,Prior,Mcmc)
{
#
# revision history:
#   changed 12/17/04 by rossi to fix bug in drawdelta when there is zero/one unit
#   in a mixture component
#   adapted to linear model by Vicky Chen 6/06
#   put in classes 3/07
#   changed a check 9/08
#
# purpose: run hierarchical linear model with mixture of normals 
#
# Arguments:
#   Data contains a list of (regdata, and possibly Z)
#      regdata is a list of lists (one list per unit)
#          regdata[[i]]=list(y,X)
#             y is a vector of observations
#             X is a length(y) x nvar matrix of values of
#               X vars including intercepts
#             Z is an nreg x nz matrix of values of variables
#               note: Z should NOT contain an intercept
#   Prior contains a list of (nu.e,ssq,deltabar,Ad,mubar,Amu,nu,V,ncomp,a) 
#      ncomp is the number of components in normal mixture
#           if elements of Prior (other than ncomp) do not exist, defaults are used
#   Mcmc contains a list of (s,c,R,keep)
#
# Output:  as list containing
#   taodraw is R/keep x nreg  array of error variances for each regression
#   Deltadraw R/keep  x nz*nvar matrix of draws of Delta, first row is initial value
#   betadraw is nreg x nvar x R/keep array of draws of betas
#   probdraw is R/keep x ncomp matrix of draws of probs of mixture components
#   compdraw is a list of list of lists (length R/keep)
#      compdraw[[rep]] is the repth draw of components for mixtures
#
# Priors:
#    tau_i ~ nu.e*ssq_i/chisq(nu.e)  tau_i is the variance of epsilon_i
#    beta_i = delta %*% z[i,] + u_i
#       u_i ~ N(mu_ind[i],Sigma_ind[i])
#       ind[i] ~multinomial(p)
#       p ~ dirichlet (a)
#           a: Dirichlet parameters for prior on p
#       delta is a k x nz array
#          delta= vec(D) ~ N(deltabar,A_d^-1)
#    mu_j ~ N(mubar,A_mu^-1(x)Sigma_j)
#    Sigma_j ~ IW(nu,V^-1)
#    ncomp is number of components
#
# MCMC parameters
#   R is number of draws
#   keep is thinning parameter, keep every keepth draw
#
#  check arguments
#
#--------------------------------------------------------------------------------------------------
#
#  create functions needed
#
append=function(l) { l=c(l,list(XpX=crossprod(l$X),Xpy=crossprod(l$X,l$y)))}
#
getvar=function(l) { 
     v=var(l$y)
     if(is.na(v)) return(1)
     if(v>0) return (v) else return (1)}
#
runiregG=
function(y,X,XpX,Xpy,sigmasq,rooti,betabar,nu,ssq){
# 
# Purpose:
#   perform one Gibbs iteration for Univ Regression Model
#   only does one iteration so can be used in both rhierLinearMixture & rhierLinearModel
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
A=crossprod(rooti)
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
#
drawDelta=
function(x,y,z,comps,deltabar,Ad){
# Z,oldbetas,ind,oldcomp,deltabar,Ad
# delta = vec(D)
#  given z and comps (z[i] gives component indicator for the ith observation, 
#   comps is a list of mu and rooti)
# y is betas: nreg x nvar
# x is Z: nreg x nz
# y = xD' + U , rows of U are indep with covs Sigma_i given by z and comps
nvar=ncol(y) #p
nz=ncol(x)   #k
xtx = matrix(0.0,nz*nvar,nz*nvar)
xty = matrix(0.0,nvar,nz) #this is the unvecced version, have to vec after sum
for(i in 1:length(comps)) {
   nobs=sum(z==i)
   if(nobs > 0) {
      if(nobs == 1) 
        { yi = matrix(y[z==i,],ncol=nvar); xi = matrix(x[z==i,],ncol=nz)}
      else
        { yi = y[z==i,]; xi = x[z==i,]}
          
      yi = t(t(yi)-comps[[i]][[1]])
      sigi = crossprod(t(comps[[i]][[2]]))
      xtx = xtx + crossprod(xi) %x% sigi
      xty = xty + (sigi %*% crossprod(yi,xi))
      }
}
xty = matrix(xty,ncol=1)

# then vec(t(D)) ~ N(V^{-1}(xty + Ad*deltabar),V^{-1}) V = (xtx+Ad)
cov=chol2inv(chol(xtx+Ad))
return(cov%*%(xty+Ad%*%deltabar) + t(chol(cov))%*%rnorm(length(deltabar)))
}
#-------------------------------------------------------------------------------------------------------
pandterm=function(message) { stop(message,call.=FALSE) }
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata, and (possibly) Z")}
  if(is.null(Data$regdata)) {pandterm("Requires Data element regdata (list of data for each unit)")}
  regdata=Data$regdata
  nreg=length(regdata)
  drawdelta=TRUE
if(is.null(Data$Z)) { cat("Z not specified",fill=TRUE); fsh() ; drawdelta=FALSE}
  else {if (nrow(Data$Z) != nreg) {pandterm(paste("Nrow(Z) ",nrow(Z),"ne number regressions ",nreg))}
      else {Z=Data$Z}}
  if(drawdelta) {
     nz=ncol(Z)
     colmeans=apply(Z,2,mean)
     if(sum(colmeans) > .00001) 
       {pandterm(paste("Z does not appear to be de-meaned: colmeans= ",colmeans))}
  }
#
# check regdata for validity
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
# check on prior
#
if(missing(Prior)) 
{pandterm("Requires Prior list argument (at least ncomp)")} 
if(is.null(Prior$nu.e)) {nu.e=3} 
   else {nu.e=Prior$nu.e}
if(is.null(Prior$ssq)) {ssq=sapply(regdata,getvar)} 
   else {ssq=Prior$ssq}
if(is.null(Prior$ncomp)) {pandterm("Requires Prior element ncomp (num of mixture components)")} else {ncomp=Prior$ncomp}
if(is.null(Prior$mubar)) {mubar=matrix(rep(0,nvar),nrow=1)} else { mubar=matrix(Prior$mubar,nrow=1)}
  if(ncol(mubar) != nvar) {pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ",ncol(mubar)))}
if(is.null(Prior$Amu)) {Amu=matrix(.01,ncol=1)} else {Amu=matrix(Prior$Amu,ncol=1)}
  if(ncol(Amu) != 1 | nrow(Amu) != 1) {pandterm("Am must be a 1 x 1 array")}
if(is.null(Prior$nu)) {nu=nvar+3}  else {nu=Prior$nu}
  if(nu < 1) {pandterm("invalid nu value")}
if(is.null(Prior$V)) {V=nu*diag(nvar)} else {V=Prior$V}
  if(sum(dim(V)==c(nvar,nvar)) !=2) pandterm("Invalid V in prior")
if(is.null(Prior$Ad) & drawdelta) {Ad=.01*diag(nvar*nz)} else {Ad=Prior$Ad}
if(drawdelta) {if(ncol(Ad) != nvar*nz | nrow(Ad) != nvar*nz) {pandterm("Ad must be nvar*nz x nvar*nz")}}
if(is.null(Prior$deltabar)& drawdelta) {deltabar=rep(0,nz*nvar)} else {deltabar=Prior$deltabar}
  if(drawdelta) {if(length(deltabar) != nz*nvar) {pandterm("deltabar must be of length nvar*nz")}}
if(is.null(Prior$a)) { a=rep(5,ncomp)} else {a=Prior$a}
if(length(a) != ncomp) {pandterm("Requires dim(a)= ncomp (no of components)")}
bada=FALSE
   for(i in 1:ncomp) { if(a[i] < 0) bada=TRUE}
  if(bada) pandterm("invalid values in a vector")
#
# check on Mcmc
#
if(missing(Mcmc)) 
  {pandterm("Requires Mcmc list argument")}
else 
   { 
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
    if(is.null(Mcmc$R)) {pandterm("Requires R argument in Mcmc list")} else {R=Mcmc$R}
    }
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting MCMC Inference for Hierarchical Linear Model:",fill=TRUE)
cat("   Normal Mixture with",ncomp,"components for first stage prior",fill=TRUE)
cat(paste("   for ",nreg," cross-sectional units"),fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("nu.e =",nu.e,fill=TRUE)
cat("nu =",nu,fill=TRUE)
cat("V ",fill=TRUE)
print(V)
cat("mubar ",fill=TRUE)
print(mubar)
cat("Amu ", fill=TRUE)
print(Amu)
cat("a ",fill=TRUE)
print(a)
if(drawdelta) 
{
   cat("deltabar",fill=TRUE)
   print(deltabar)
   cat("Ad",fill=TRUE)
   print(Ad)
}
cat(" ",fill=TRUE)
cat("MCMC Parms: ",fill=TRUE)
cat(paste(" R= ",R," keep= ",keep),fill=TRUE)
cat("",fill=TRUE)
#
# allocate space for draws
#
taudraw=matrix(double(floor(R/keep)*nreg),ncol=nreg)
if(drawdelta) Deltadraw=matrix(double((floor(R/keep))*nz*nvar),ncol=nz*nvar)
betadraw=array(double((floor(R/keep))*nreg*nvar),dim=c(nreg,nvar,floor(R/keep)))
probdraw=matrix(double((floor(R/keep))*ncomp),ncol=ncomp)
oldbetas=matrix(double(nreg*nvar),ncol=nvar)
oldcomp=NULL
compdraw=NULL
#
#  initialize values
#
#  Create XpX elements of regdata and initialize tau
#
regdata=lapply(regdata,append)
tau=sapply(regdata,getvar)
#
# set initial values for the indicators
#     ind is of length(nreg) and indicates which mixture component this obs
#     belongs to.
#
ind=NULL
ninc=floor(nreg/ncomp)
for (i in 1:(ncomp-1)) {ind=c(ind,rep(i,ninc))}
if(ncomp != 1) {ind = c(ind,rep(ncomp,nreg-length(ind)))} else {ind=rep(1,nreg)}
#
# initialize delta
#
if (drawdelta) olddelta=rep(0,nz*nvar)
#
# initialize probs
#
oldprob=rep(1/ncomp,ncomp)
#
# initialize comps
#
tcomp=list(list(mu=rep(0,nvar),rooti=diag(nvar)))
oldcomp=rep(tcomp,ncomp)
#
#	start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
fsh()
for(rep in 1:R)
{
   # first draw comps,ind,p | {beta_i}, delta
   #        ind,p need initialization comps is drawn first in sub-Gibbs
   if(drawdelta) 
      {mgout=rmixGibbs(oldbetas-Z%*%t(matrix(olddelta,ncol=nz)),
      mubar,Amu,nu,V,a,oldprob,ind,oldcomp)}
   else
      {mgout=rmixGibbs(oldbetas,
      mubar,Amu,nu,V,a,oldprob,ind,oldcomp)}
   oldprob=mgout[[1]]
   oldcomp=mgout[[3]]
   ind=mgout[[2]]
   # now draw delta | {beta_i}, ind, comps
   if(drawdelta) {olddelta=drawDelta(Z,oldbetas,ind,oldcomp,deltabar,Ad)}
   #
   #  loop over all regression equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
   #
      for (reg in 1:nreg) 
      {
         rootpi=oldcomp[[ind[reg]]]$rooti
         #  note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
         if(drawdelta) {
            betabar=oldcomp[[ind[reg]]]$mu+matrix(olddelta,ncol=nz)%*%as.vector(Z[reg,])}
         else {
            betabar=oldcomp[[ind[reg]]]$mu }
      regout=runiregG(regdata[[reg]]$y,regdata[[reg]]$X,regdata[[reg]]$XpX,
                regdata[[reg]]$Xpy,tau[reg],rootpi,betabar,nu.e,ssq[reg])
      oldbetas[reg,]=regout$betadraw
      tau[reg]=regout$sigmasqdraw
      }
   #
   #       print time to completion and draw # every 100th draw
   #
   if(((rep/100)*100) ==(floor(rep/100)*100))
     {ctime=proc.time()[3]
      timetoend=((ctime-itime)/rep)*(R+1-rep)
      cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
      fsh()}
   #
   #       save every keepth draw
   #
   mkeep=rep/keep
   if((mkeep*keep) == (floor(mkeep)*keep))
      { taudraw[mkeep,]=tau
        betadraw[,,mkeep]=oldbetas 
        probdraw[mkeep,]=oldprob
        if(drawdelta) Deltadraw[mkeep,]=olddelta
        compdraw[[mkeep]]=oldcomp }
        
}
ctime=proc.time()[3]
cat(" Total Time Elapsed: ",round((ctime-itime)/60,2),fill=TRUE)
attributes(taudraw)$class=c("bayesm.mat","mcmc")
attributes(taudraw)$mcpar=c(1,R,keep)
if(drawdelta){
   attributes(Deltadraw)$class=c("bayesm.mat","mcmc")
   attributes(Deltadraw)$mcpar=c(1,R,keep)}
attributes(betadraw)$class=c("bayesm.hcoef")
nmix=list(probdraw=probdraw,zdraw=NULL,compdraw=compdraw)
attributes(nmix)$class="bayesm.nmix"
if(drawdelta) 
   {return(list(taudraw=taudraw,Deltadraw=Deltadraw,betadraw=betadraw,nmix=nmix))} 
else 
   {return(list(taudraw=taudraw,betadraw=betadraw,nmix=nmix))}
}
