rhierMnlRwMixture=
function(Data,Prior,Mcmc)
{
#
# revision history:
#   changed 12/17/04 by rossi to fix bug in drawdelta when there is zero/one unit
#   in a mixture component
#   added loglike output, changed to reflect new argument order in llmnl, mnlHess 9/05
#   changed weighting scheme to (1-w)logl_i + w*Lbar (normalized) 12/05
#   3/07 added classes
#
# purpose: run hierarchical mnl logit model with mixture of normals 
#   using RW and cov(RW inc) = (hess_i + Vbeta^-1)^-1
#   uses normal approximation to pooled likelihood
#
# Arguments:
#   Data contains a list of (p,lgtdata, and possibly Z)
#      p is number of choice alternatives
#      lgtdata is a list of lists (one list per unit)
#          lgtdata[[i]]=list(y,X)
#             y is a vector indicating alternative chosen
#               integers 1:p indicate alternative
#             X is a length(y)*p x nvar matrix of values of
#               X vars including intercepts
#             Z is an length(lgtdata) x nz matrix of values of variables
#               note: Z should NOT contain an intercept
#   Prior contains a list of (deltabar,Ad,mubar,Amu,nu,V,ncomp) 
#      ncomp is the number of components in normal mixture
#           if elements of Prior (other than ncomp) do not exist, defaults are used
#   Mcmc contains a list of (s,c,R,keep)
#
# Output:  as list containing
#   Deltadraw R/keep  x nz*nvar matrix of draws of Delta, first row is initial value
#   betadraw is nlgt x nvar x R/keep array of draws of betas
#   probdraw is R/keep x ncomp matrix of draws of probs of mixture components
#   compdraw is a list of list of lists (length R/keep)
#      compdraw[[rep]] is the repth draw of components for mixtures
#   loglike  log-likelikelhood at each kept draw
#
# Priors:
#    beta_i = D %*% z[i,] + u_i
#       u_i ~ N(mu_ind[i],Sigma_ind[i])
#       ind[i] ~multinomial(p)
#       p ~ dirichlet (a)
#       D is a k x nz array
#          delta= vec(D) ~ N(deltabar,A_d^-1)
#    mu_j ~ N(mubar,A_mu^-1(x)Sigma_j)
#    Sigma_j ~ IW(nu,V^-1)
#    ncomp is number of components
#
# MCMC parameters
#   s is the scaling parameter for the RW inc covariance matrix; s^2 Var is inc cov
#      matrix
#   w is parameter for weighting function in fractional likelihood
#      w is the weight on the normalized pooled likelihood 
#   R is number of draws
#   keep is thinning parameter, keep every keepth draw
#
#  check arguments
#
pandterm=function(message) { stop(message,call.=FALSE) }
if(missing(Data)) {pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")}
  if(is.null(Data$p)) {pandterm("Requires Data element p (# chce alternatives)") }
  p=Data$p
  if(is.null(Data$lgtdata)) {pandterm("Requires Data element lgtdata (list of data for each unit)")}
  lgtdata=Data$lgtdata
  nlgt=length(lgtdata)
  drawdelta=TRUE
if(is.null(Data$Z)) { cat("Z not specified",fill=TRUE); fsh() ; drawdelta=FALSE}
  else {if (nrow(Data$Z) != nlgt) {pandterm(paste("Nrow(Z) ",nrow(Z),"ne number logits ",nlgt))}
      else {Z=Data$Z}}
  if(drawdelta) {
     nz=ncol(Z)
     colmeans=apply(Z,2,mean)
     if(sum(colmeans) > .00001) 
       {pandterm(paste("Z does not appear to be de-meaned: colmeans= ",colmeans))}
  }
  
#
# check lgtdata for validity
#
ypooled=NULL
Xpooled=NULL
if(!is.null(lgtdata[[1]]$X)) {oldncol=ncol(lgtdata[[1]]$X)}
for (i in 1:nlgt) 
{
    if(is.null(lgtdata[[i]]$y)) {pandterm(paste("Requires element y of lgtdata[[",i,"]]"))}
    if(is.null(lgtdata[[i]]$X)) {pandterm(paste("Requires element X of lgtdata[[",i,"]]"))}
    ypooled=c(ypooled,lgtdata[[i]]$y)
    nrowX=nrow(lgtdata[[i]]$X)
    if((nrowX/p) !=length(lgtdata[[i]]$y)) {pandterm(paste("nrow(X) ne p*length(yi); exception at unit",i))}
    newncol=ncol(lgtdata[[i]]$X)
    if(newncol != oldncol) {pandterm(paste("All X elements must have same # of cols; exception at unit",i))}
    Xpooled=rbind(Xpooled,lgtdata[[i]]$X)
    oldncol=newncol
}
nvar=ncol(Xpooled)
levely=as.numeric(levels(as.factor(ypooled)))
if(length(levely) != p) {pandterm(paste("y takes on ",length(levely)," values -- must be = p"))}
bady=FALSE
for (i in 1:p )
{
    if(levely[i] != i) bady=TRUE
}
cat("Table of Y values pooled over all units",fill=TRUE)
print(table(ypooled))
if (bady) 
  {pandterm("Invalid Y")}
#
# check on prior
#
if(missing(Prior)) 
{pandterm("Requires Prior list argument (at least ncomp)")} 
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
   for(i in 1:ncomp) { if(a[i] < 1) bada=TRUE}
  if(bada) pandterm("invalid values in a vector")
#
# check on Mcmc
#
if(missing(Mcmc)) 
  {pandterm("Requires Mcmc list argument")}
else 
   { 
    if(is.null(Mcmc$s)) {s=2.93/sqrt(nvar)} else {s=Mcmc$s}
    if(is.null(Mcmc$w)) {w=.1}  else {w=Mcmc$w}
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
    if(is.null(Mcmc$R)) {pandterm("Requires R argument in Mcmc list")} else {R=Mcmc$R}
    }
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting MCMC Inference for Hierarchical Logit:",fill=TRUE)
cat("   Normal Mixture with",ncomp,"components for first stage prior",fill=TRUE)
cat(paste("  ",p," alternatives; ",nvar," variables in X"),fill=TRUE)
cat(paste("   for ",nlgt," cross-sectional units"),fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
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
cat(paste("s=",round(s,3)," w= ",w," R= ",R," keep= ",keep),fill=TRUE)
cat("",fill=TRUE)
#
# allocate space for draws
#
if(drawdelta) Deltadraw=matrix(double((floor(R/keep))*nz*nvar),ncol=nz*nvar)
betadraw=array(double((floor(R/keep))*nlgt*nvar),dim=c(nlgt,nvar,floor(R/keep)))
probdraw=matrix(double((floor(R/keep))*ncomp),ncol=ncomp)
oldbetas=matrix(double(nlgt*nvar),ncol=nvar)
oldll=double(nlgt)
loglike=double(floor(R/keep))
oldcomp=NULL
compdraw=NULL
#--------------------------------------------------------------------------------------------------
#
#  create functions needed
#
llmnlFract=
function(beta,y,X,betapooled,rootH,w,wgt){
z=as.vector(rootH%*%(beta-betapooled))
return((1-w)*llmnl(beta,y,X)+w*wgt*(-.5*(z%*%z)))
}

mnlRwMetropOnce=
function(y,X,oldbeta,oldll,s,inc.root,betabar,rootpi){ 
#
# function to execute rw metropolis for the MNL
# y is n vector with element = 1,...,j indicating which alt chosen
# X is nj x k matrix of xvalues for each of j alt on each of n occasions
# RW increments are N(0,s^2*t(inc.root)%*%inc.root)
# prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
#	inc.root, rootpi are upper triangular
#	this means that we are using the UL decomp of Sigma^-1 for prior 
# oldbeta is the current
     stay=0
     betac=oldbeta + s*t(inc.root)%*%(matrix(rnorm(ncol(X)),ncol=1))
     cll=llmnl(betac,y,X)
     clpost=cll+lndMvn(betac,betabar,rootpi)
     ldiff=clpost-oldll-lndMvn(oldbeta,betabar,rootpi)
     alpha=min(1,exp(ldiff))
     if(alpha < 1) {unif=runif(1)} else {unif=0}
     if (unif <= alpha)
             {betadraw=betac; oldll=cll}
           else
             {betadraw=oldbeta; stay=1}
return(list(betadraw=betadraw,stay=stay,oldll=oldll))
}
drawDelta=
function(x,y,z,comps,deltabar,Ad){
# delta = vec(D)
#  given z and comps (z[i] gives component indicator for the ith observation, 
#   comps is a list of mu and rooti)
#y is n x p
#x is n x k
#y = xD' + U , rows of U are indep with covs Sigma_i given by z and comps
p=ncol(y)
k=ncol(x)
xtx = matrix(0.0,k*p,k*p)
xty = matrix(0.0,p,k) #this is the unvecced version, have to vec after sum
for(i in 1:length(comps)) {
   nobs=sum(z==i)
   if(nobs > 0) {
      if(nobs == 1) 
        { yi = matrix(y[z==i,],ncol=p); xi = matrix(x[z==i,],ncol=k)}
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
#
# intialize compute quantities for Metropolis
#
cat("initializing Metropolis candidate densities for ",nlgt," units ...",fill=TRUE)
fsh()
#
#  now go thru and computed fraction likelihood estimates and hessians
#
#       Lbar=log(pooled likelihood^(n_i/N))
#
#       fraction loglike = (1-w)*loglike_i + w*Lbar
#
betainit=c(rep(0,nvar))
#
#  compute pooled optimum
#
out=optim(betainit,llmnl,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-6), 
     X=Xpooled,y=ypooled)
betapooled=out$par
H=mnlHess(betapooled,ypooled,Xpooled)
rootH=chol(H)
for (i in 1:nlgt) 
{
   wgt=length(lgtdata[[i]]$y)/length(ypooled)
   out=optim(betapooled,llmnlFract,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-4), 
   X=lgtdata[[i]]$X,y=lgtdata[[i]]$y,betapooled=betapooled,rootH=rootH,w=w,wgt=wgt)
   if(out$convergence == 0) 
     { hess=mnlHess(out$par,lgtdata[[i]]$y,lgtdata[[i]]$X)
       lgtdata[[i]]=c(lgtdata[[i]],list(converge=1,betafmle=out$par,hess=hess)) }
   else
     { lgtdata[[i]]=c(lgtdata[[i]],list(converge=0,betafmle=c(rep(0,nvar)),
        hess=diag(nvar))) }
   oldbetas[i,]=lgtdata[[i]]$betafmle
   if(i%%50 ==0) cat("  completed unit #",i,fill=TRUE)
   fsh()
}
#
#  initialize values
#
# set initial values for the indicators
#     ind is of length(nlgt) and indicates which mixture component this obs
#     belongs to.
#
ind=NULL
ninc=floor(nlgt/ncomp)
for (i in 1:(ncomp-1)) {ind=c(ind,rep(i,ninc))}
if(ncomp != 1) {ind = c(ind,rep(ncomp,nlgt-length(ind)))} else {ind=rep(1,nlgt)}
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
   #  loop over all lgt equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
   #
      for (lgt in 1:nlgt) 
      {
         rootpi=oldcomp[[ind[lgt]]]$rooti
         #  note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
         if(drawdelta) {
            betabar=oldcomp[[ind[lgt]]]$mu+matrix(olddelta,ncol=nz)%*%as.vector(Z[lgt,])}
         else {
            betabar=oldcomp[[ind[lgt]]]$mu }
         if (rep == 1) 
            { oldll[lgt]=llmnl(oldbetas[lgt,],lgtdata[[lgt]]$y,lgtdata[[lgt]]$X)}  
         #   compute inc.root
         inc.root=chol(chol2inv(chol(lgtdata[[lgt]]$hess+rootpi%*%t(rootpi))))
         metropout=mnlRwMetropOnce(lgtdata[[lgt]]$y,lgtdata[[lgt]]$X,oldbetas[lgt,],
                                   oldll[lgt],s,inc.root,betabar,rootpi)      
         oldbetas[lgt,]=metropout$betadraw
         oldll[lgt]=metropout$oldll
      }
   #
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
      { betadraw[,,mkeep]=oldbetas 
        probdraw[mkeep,]=oldprob
        loglike[mkeep]=sum(oldll)
        if(drawdelta) Deltadraw[mkeep,]=olddelta
        compdraw[[mkeep]]=oldcomp }
        
}
ctime=proc.time()[3]
cat(" Total Time Elapsed: ",round((ctime-itime)/60,2),fill=TRUE)
if(drawdelta){
   attributes(Deltadraw)$class=c("bayesm.mat","mcmc")
   attributes(Deltadraw)$mcpar=c(1,R,keep)}
attributes(betadraw)$class=c("bayesm.hcoef")
nmix=list(probdraw=probdraw,zdraw=NULL,compdraw=compdraw)
attributes(nmix)$class="bayesm.nmix"
if(drawdelta) 
   {return(list(Deltadraw=Deltadraw,betadraw=betadraw,nmix=nmix,loglike=loglike))} 
else 
   {return(list(betadraw=betadraw,nmix=nmix,loglike=loglike))}
}
