rhierMnlDP=
function(Data,Prior,Mcmc)
{
#
#  created 3/08 by Rossi from rhierMnlRwMixture adding DP draw for to replace finite mixture of normals
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
#   Prior contains a list of (deltabar,Ad,lambda_hyper,Prioralpha)
#       alpha: starting value
#       lambda_hyper: hyperparms of prior on lambda
#       Prioralpha: hyperparms of alpha prior; a list of (Istarmin,Istarmax,power)
#       if elements of the prior don't exist, defaults are assumed
#   Mcmc contains a list of (s,c,R,keep)
#
# Output:  as list containing
#   Deltadraw R/keep  x nz*nvar matrix of draws of Delta, first row is initial value
#   betadraw is nlgt x nvar x R/keep array of draws of betas
#   probdraw is R/keep x 1 matrix of draws of probs of mixture components
#   compdraw is a list of list of lists (length R/keep)
#      compdraw[[rep]] is the repth draw of components for mixtures
#   loglike  log-likelikelhood at each kept draw
#
# Priors:
#    beta_i = D %*% z[i,] + u_i
#       vec(D)~N(deltabar)
#       u_i ~ N(theta_i)
#       theta_i~G
#       G|lambda,alpha ~ DP(G|G0(lambda),alpha)
#
#        lambda:
#           G0 ~ N(mubar,Sigma (x) Amu^-1)
#           mubar=vec(mubar)
#           Sigma ~ IW(nu,nu*v*I)  note: mode(Sigma)=nu/(nu+2)*v*I
#           mubar=0
#           amu is uniform on grid specified by alim
#           nu is log uniform, nu=d-1+exp(Z) z is uniform on seq defined bvy nulim
#           v is uniform on sequence specificd by vlim
#
#        Prioralpha:
#           alpha ~ (1-(alpha-alphamin)/(alphamax-alphamin))^power
#           alphamin=exp(digamma(Istarmin)-log(gamma+log(N)))
#           alphamax=exp(digamma(Istarmax)-log(gamma+log(N)))
#           gamma= .5772156649015328606
#
# MCMC parameters
#   s is the scaling parameter for the RW inc covariance matrix; s^2 Var is inc cov
#      matrix
#   w is parameter for weighting function in fractional likelihood
#      w is the weight on the normalized pooled likelihood 
#   R is number of draws
#   keep is thinning parameter, keep every keepth draw
#--------------------------------------------------------------------------------------------------
#
#  create functions needed
#
rDPGibbs1= 
function(y,theta,thetaStar,indic,lambda,alpha,Prioralpha,lambda_hyper,maxuniq,gridsize){
#
#  revision history:
#   created from rDPGibbs by Rossi 3/08
#
#  do one draw of DP Gibbs sampler with norma base
#
# Model:
#        y_i ~ N(y|thetai)
#        thetai|G ~ G
#        G|lambda,alpha ~ DP(G|G0(lambda),alpha)
#
# Priors:
#        alpha: starting value
#
#        lambda:
#           G0 ~ N(mubar,Sigma (x) Amu^-1)
#           mubar=vec(mubar)
#           Sigma ~ IW(nu,nu*V) V=v*I  note: mode(Sigma)=nu/(nu+2)*v*I
#           mubar=0
#           amu is uniform on grid specified by alim
#           nu is log uniform, nu=d-1+exp(Z) z is uniform on seq defined bvy nulim
#           v is uniform on sequence specificd by vlim
#
#        Prioralpha:
#           alpha ~ (1-(alpha-alphamin)/(alphamax-alphamin))^power
#           alphamin=exp(digamma(Istarmin)-log(gamma+log(N)))
#           alphamax=exp(digamma(Istarmax)-log(gamma+log(N)))
#           gamma= .5772156649015328606
#
#
# output:
#   theta - list of thetas for each "obs"
#   ind - vector of indicators for which observations are associated with which comp in thetaStar
#   thetaStar - list of unique normal component parms
#   lambda  - list of of (a,nu,V)
#   alpha 
#   thetaNp1 - one draw from predictive given thetaStar, lambda,alphama
#
# define needed functions
#
# -----------------------------------------------------------------------------------------
#
q0=
   function(y,lambda,eta){
#
# function to compute a vector of int f(y[i]|theta) p(theta|lambda)dlambda
#     here p(theta|lambda) is G0 the base prior
#
# implemented for a multivariate normal data density and standard conjugate
# prior:
#    theta=list(mu,Sigma)
#    f(y|theta,eta) is N(mu,Sigma)
#    lambda=list(mubar,Amu,nu,V)
#       mu|Sigma ~ N(mubar,Sigma (x) Amu^-1)
#       Sigma ~ IW(nu,V)
#
# arguments:
#    Y is n x k matrix of observations
#    lambda=list(mubar,Amu,nu,V)
#    eta is not used
# 
# output:
#    vector of q0 values for each obs (row of Y)
#
# p. rossi 12/05
#
# here y is matrix of observations (each row is an obs)

mubar=lambda$mubar; nu=lambda$nu ; Amu=lambda$Amu; V=lambda$V
k=ncol(y)
R=chol(V)
logdetR=sum(log(diag(R)))
if (k > 1) 
  {lnk1k2=(k/2)*log(2)+log((nu-k)/2)+lgamma((nu-k)/2)-lgamma(nu/2)+sum(log(nu/2-(1:(k-1))/2))}
else
  {lnk1k2=(k/2)*log(2)+log((nu-k)/2)+lgamma((nu-k)/2)-lgamma(nu/2)}
constant=-(k/2)*log(2*pi)+(k/2)*log(Amu/(1+Amu)) + lnk1k2 + nu*logdetR
#
# note: here we are using the fact that |V + S_i | = |R|^2 (1 + v_i'v_i)
#       where v_i = sqrt(Amu/(1+Amu))*t(R^-1)*(y_i-mubar), R is chol(V)
#
#       and S_i = Amu/(1+Amu) * (y_i-mubar)(y_i-mubar)'
#
mat=sqrt(Amu/(1+Amu))*t(backsolve(R,diag(ncol(y))))%*%(t(y)-mubar)
vivi=colSums(mat^2)

lnq0v=constant-((nu+1)/2)*(2*logdetR+log(1+vivi))

return(exp(lnq0v))
}
# ----------------------------------------------------------------------------------------------
   rmultinomF=function(p) {
       return(sum(runif(1) > cumsum(p))+1)
   }
# -----------------------------------------------------------------------------------------------
alphaD=
   function(Prioralpha,Istar,gridsize){
#
#  function to draw alpha using prior, p(alpha)= (1-(alpha-alphamin)/(alphamax-alphamin))**power
#
   power=Prioralpha$power
   alphamin=Prioralpha$alphamin
   alphamax=Prioralpha$alphamax
   n=Prioralpha$n
   alpha=seq(from=alphamin,to=(alphamax-0.000001),len=gridsize)
   lnprob=Istar*log(alpha) + lgamma(alpha) - lgamma(n+alpha) + 
          power*log(1-(alpha-alphamin)/(alphamax-alphamin))
   lnprob=lnprob-median(lnprob)
   probs=exp(lnprob)
   probs=probs/sum(probs)
   return(alpha[rmultinomF(probs)])
}  


#
# ------------------------------------------------------------------------------------------
#
yden=
   function(thetaStar,y,eta){
#
# function to compute f(y | theta) 
# computes f for all values of theta in theta list of lists
#
# arguments:
#   thetaStar is a list of lists.  thetaStar[[i]] is a list with components, mu, rooti
#   y |theta[[i]] ~ N(mu,(rooti %*% t(rooti))^-1)  rooti is inverse of Chol root of Sigma
#   eta is not used
#
# output:
#   length(thetaStar) x n array of values of f(y[j,]|thetaStar[[i]]
# 

nunique=length(thetaStar)
n=nrow(y)
ydenmat=matrix(double(n*nunique),ncol=n)
k=ncol(y)
for(i in 1:nunique){

   # now compute vectorized version of lndMvn 
   # compute y_i'RIRI'y_i for all i
   #
   mu=thetaStar[[i]]$mu; rooti=thetaStar[[i]]$rooti
   quads=colSums((crossprod(rooti,(t(y)-mu)))^2)
   ydenmat[i,]=exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads)
   
}
return(ydenmat)
}

#
# -----------------------------------------------------------------------------------------
#
GD=
   function(lambda){
#
# function to draw from prior for Multivariate Normal Model
#
# mu|Sigma ~ N(mubar,Sigma x Amu^-1)
# Sigma ~ IW(nu,V)
#
# note: we must insure that mu is a vector to use most efficient
#       lndMvn routine
#
nu=lambda$nu
V=lambda$V
mubar=lambda$mubar
Amu=lambda$Amu
k=length(mubar)
Sigma=rwishart(nu,chol2inv(chol(lambda$V)))$IW
root=chol(Sigma)
mu=mubar+(1/sqrt(Amu))*t(root)%*%matrix(rnorm(k),ncol=1)
return(list(mu=as.vector(mu),rooti=backsolve(root,diag(k))))
}

#
# -------------------------------------------------------------------------------------------
#
thetaD=
   function(y,lambda,eta){
#
# function to draw from posterior of theta given data y and base prior G0(lambda)
#
# here y ~ N(mu,Sigma)
# theta = list(mu=mu,rooti=chol(Sigma)^-1)
# mu|Sigma ~ N(mubar,Sigma (x) Amu-1)
# Sigma ~ IW(nu,V)
#
# arguments: 
#   y is n x k matrix of obs
#   lambda is list(mubar,Amu,nu,V)
#   eta is not used
# output:
#   one draw of theta, list(mu,rooti)
#        Sigma=inv(rooti)%*%t(inv(rooti))
#
# note: we assume that y is a matrix. if there is only one obs, y is a 1 x k matrix
#
rout=rmultireg(y,matrix(c(rep(1,nrow(y))),ncol=1),matrix(lambda$mubar,nrow=1),matrix(lambda$Amu,ncol=1),
       lambda$nu,lambda$V)
return(list(mu=as.vector(rout$B),rooti=backsolve(chol(rout$Sigma),diag(ncol(y)))))
}

#
# --------------------------------------------------------------------------------------------
# load a faster version of lndMvn
# note: version of lndMvn below assumes x,mu is a vector!
lndMvn=function (x, mu, rooti){
    return(-(length(x)/2) * log(2 * pi) - 0.5 * sum(((x-mu)%*%rooti)**2) + sum(log(diag(rooti))))
}
# -----------------------------------------------------------------------------------------
   lambdaD=function(lambda,thetaStar,alim=c(.01,2),nulim=c(.01,2),vlim=c(.1,5),gridsize=20){
#
# revision history
#  p. rossi 7/06
#  vectorized 1/07
#  changed 2/08 to paramaterize V matrix of IW prior to nu*v*I; then mode of Sigma=nu/(nu+2)vI
#      this means that we have a reparameterization to v* = nu*v
#
#  function to draw (nu, v, a) using uniform priors
#
#  theta_j=(mu_j,Sigma_j)  mu_j~N(0,Sigma_j/a)  Sigma_j~IW(nu,vI)
#           recall E[Sigma]= vI/(nu-dim-1)
#
# define functions needed
# ----------------------------------------------------------------------------------------------
   rmultinomF=function(p) {
       return(sum(runif(1) > cumsum(p))+1)
   }
echo=function(lst){return(t(lst[[2]]))}
rootiz=function(lst){crossprod(lst[[2]],lst[[1]])}
#
# ------------------------------------------------------------------------------------------

   d=length(thetaStar[[1]]$mu)
   Istar=length(thetaStar)
   aseq=seq(from=alim[1],to=alim[2],len=gridsize)
   nuseq=d-1+exp(seq(from=nulim[1],to=nulim[2],len=gridsize)) # log uniform grid
   vseq=seq(from=vlim[1],to=vlim[2],len=gridsize)
#
# extract needed info from thetaStar list
#
   out=double(Istar*d*d)
   out=sapply(thetaStar,echo)
   dim(out)=c(d,Istar*d) # out has the rootis in form: [t(rooti_1), t(rooti_2), ...,t(rooti_Istar)]
   sumdiagriri=sum(colSums(out^2)) #  sum_j tr(rooti_j%*%t(rooti_j))
#   now get diagonals of rooti
   ind=cbind(c(1:(d*Istar)),rep((1:d),Istar))
   out=t(out)
   sumlogdiag=sum(log(out[ind]))
   rimu=sapply(thetaStar,rootiz) # columns of rimu contain t(rooti_j)%*%mu_j
   dim(rimu)=c(d,Istar)
   sumquads=sum(colSums(rimu^2)) 
#  
#  draw a  (conditionally indep of nu,v given theta_j)
    lnprob=double(length(aseq))
    lnprob=Istar*(-(d/2)*log(2*pi))-.5*aseq*sumquads+Istar*d*log(sqrt(aseq))+sumlogdiag
    lnprob=lnprob-max(lnprob)+200
    probs=exp(lnprob)
    probs=probs/sum(probs)
    adraw=aseq[rmultinomF(probs)]
#
#   draw nu given v
#
    V=lambda$V
    lnprob=double(length(nuseq))
    arg=rep(c(1:d),gridsize)
    dim(arg)=c(d,gridsize)
    arg=t(arg)
    arg=(nuseq+1-arg)/2
    lnprob=-Istar*log(2)*d/2*nuseq - Istar*rowSums(lgamma(arg)) + 
            Istar*d*log(sqrt(V[1,1]))*nuseq + sumlogdiag*nuseq
    lnprob=lnprob-max(lnprob)+200
    probs=exp(lnprob)
    probs=probs/sum(probs)
    nudraw=nuseq[rmultinomF(probs)]
#
#   draw v given nu 
#
    lnprob=double(length(vseq))
    lnprob=Istar*nudraw*d*log(sqrt(vseq*nudraw))-.5*sumdiagriri*vseq*nudraw
    lnprob=lnprob-max(lnprob)+200
    probs=exp(lnprob)
    probs=probs/sum(probs)
    vdraw=vseq[rmultinomF(probs)]
#
#   put back into lambda
#
    return(list(mubar=c(rep(0,d)),Amu=adraw,nu=nudraw,V=nudraw*vdraw*diag(d)))
}
pandterm=function(message) { stop(message,call.=FALSE) }
# -----------------------------------------------------------------------------------------

for(rep in 1:1)		#note: we only do one loop!
{
   n = length(theta)

   eta=NULL    # note eta is not used
   thetaNp1=NULL
   q0v = q0(y,lambda,eta)   # now that we draw lambda we need to recompute q0v each time

   p=c(rep(1/(alpha+(n-1)),n-1),alpha/(alpha+(n-1)))

   nunique=length(thetaStar)
  
   if(nunique > maxuniq ) { pandterm("maximum number of unique thetas exceeded")} 
   ydenmat=matrix(double(maxuniq*n),ncol=n) 
   ydenmat[1:nunique,]=yden(thetaStar,y,eta)
   #  ydenmat is a length(thetaStar) x n array of density values given f(y[j,] | thetaStar[[i]]
   #  note: due to remix step (below) we must recompute ydenmat each time!

   # use .Call to draw theta list
   out= .Call("thetadraw",y,ydenmat,indic,q0v,p,theta,lambda,eta=eta,
                  thetaD=thetaD,yden=yden,maxuniq,nunique,new.env()) 

   # theta has been modified by thetadraw so we need to recreate thetaStar
   thetaStar=unique(theta)
   nunique=length(thetaStar)

   #thetaNp1 and remix
   probs=double(nunique+1)
   for(j in 1:nunique) {
       ind = which(sapply(theta,identical,thetaStar[[j]]))
       probs[j]=length(ind)/(alpha+n) 
       new_utheta=thetaD(y[ind,,drop=FALSE],lambda,eta) 
       for(i in seq(along=ind)) {theta[[ind[i]]]=new_utheta}
       indic[ind]=j
       thetaStar[[j]]=new_utheta
   }
   probs[nunique+1]=alpha/(alpha+n)
   ind=rmultinomF(probs)
   if(ind==length(probs)) {
      thetaNp1=GD(lambda)
   } else {
      thetaNp1=thetaStar[[ind]]
   }

   # draw alpha
   alpha=alphaD(Prioralpha,nunique,gridsize=gridsize)
   
   # draw lambda
   lambda=lambdaD(lambda,thetaStar,alim=lambda_hyper$alim,nulim=lambda_hyper$nulim,
             vlim=lambda_hyper$vlim,gridsize=gridsize)

}
#   note indic is the vector of indicators for each obs correspond to which thetaStar
return(list(theta=theta,thetaStar=thetaStar,thetaNp1=thetaNp1,alpha=alpha,lambda=lambda,ind=indic))
}
#--------------------------------------------------------------------------------------------------

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
alimdef=c(.01,2)
nulimdef=c(.01,3)
vlimdef=c(.1,4)
if(missing(Prior)) {Prior=NULL}

if(is.null(Prior$lambda_hyper)) {lambda_hyper=list(alim=alimdef,nulim=nulimdef,vlim=vlimdef)}
   else {lambda_hyper=Prior$lambda_hyper;
       if(is.null(lambda_hyper$alim)) {lambda_hyper$alim=alimdef}
       if(is.null(lambda_hyper$nulim)) {lambda_hyper$nulim=nulimdef} 
       if(is.null(lambda_hyper$vlim)) {lambda_hyper$vlim=vlimdef}
       }
if(is.null(Prior$Prioralpha)) {Prioralpha=list(Istarmin=1,Istarmax=min(50,0.1*nlgt),power=0.8)}
   else {Prioralpha=Prior$Prioralpha;
       if(is.null(Prioralpha$Istarmin)) {Prioralpha$Istarmin=1} else {Prioralpha$Istarmin=Prioralpha$Istarmin}
       if(is.null(Prioralpha$Istarmax)) 
       {Prioralpha$Istarmax=min(50,0.1*nlgt)} else {Prioralpha$Istarmax=Prioralpha$Istarmax}
      if(is.null(Prioralpha$power)) {Prioralpha$power=0.8}
   }	
gamma= .5772156649015328606
Prioralpha$alphamin=exp(digamma(Prioralpha$Istarmin)-log(gamma+log(nlgt)))
Prioralpha$alphamax=exp(digamma(Prioralpha$Istarmax)-log(gamma+log(nlgt)))
Prioralpha$n=nlgt
#
# check Prior arguments for valdity
#
if(lambda_hyper$alim[1]<0) {pandterm("alim[1] must be >0")}
if(lambda_hyper$nulim[1]<0) {pandterm("nulim[1] must be >0")}
if(lambda_hyper$vlim[1]<0) {pandterm("vlim[1] must be >0")}
if(Prioralpha$Istarmin <1){pandterm("Prioralpha$Istarmin must be >= 1")}
if(Prioralpha$Istarmax <= Prioralpha$Istarmin){pandterm("Prioralpha$Istarmin must be < Prioralpha$Istarmax")}	

if(is.null(Prior$Ad) & drawdelta) {Ad=.01*diag(nvar*nz)} else {Ad=Prior$Ad}
if(drawdelta) {if(ncol(Ad) != nvar*nz | nrow(Ad) != nvar*nz) {pandterm("Ad must be nvar*nz x nvar*nz")}}
if(is.null(Prior$deltabar)& drawdelta) {deltabar=rep(0,nz*nvar)} else {deltabar=Prior$deltabar}
  if(drawdelta) {if(length(deltabar) != nz*nvar) {pandterm("deltabar must be of length nvar*nz")}}
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
    if(is.null(Mcmc$maxuniq)) {maxuniq=200} else {keep=Mcmc$maxuniq}
    if(is.null(Mcmc$gridsize)) {gridsize=20} else {gridsize=Mcmc$gridsize}
    if(is.null(Mcmc$R)) {pandterm("Requires R argument in Mcmc list")} else {R=Mcmc$R}
    }
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting MCMC Inference for Hierarchical Logit:",fill=TRUE)
cat("   Dirichlet Process Prior",fill=TRUE)
cat(paste("  ",p," alternatives; ",nvar," variables in X"),fill=TRUE)
cat(paste("   for ",nlgt," cross-sectional units"),fill=TRUE)
cat(" ",fill=TRUE)
cat(" Prior Parms: ",fill=TRUE)
cat("  G0 ~ N(mubar,Sigma (x) Amu^-1)",fill=TRUE)
cat("   mubar = ",0,fill=TRUE)
cat("   Sigma ~ IW(nu,nu*v*I)",fill=TRUE)
cat("   Amu ~ uniform[",lambda_hyper$alim[1],",",lambda_hyper$alim[2],"]",fill=TRUE)
cat("   nu ~ uniform on log grid  [",nvar-1+exp(lambda_hyper$nulim[1]),
             ",",nvar-1+exp(lambda_hyper$nulim[2]),"]",fill=TRUE)
cat("   v ~ uniform[",lambda_hyper$vlim[1],",",lambda_hyper$vlim[2],"]",fill=TRUE)
cat(" ",fill=TRUE)
cat("  alpha ~ (1-(alpha-alphamin)/(alphamax-alphamin))^power",fill=TRUE)
cat("   Istarmin = ",Prioralpha$Istarmin,fill=TRUE)
cat("   Istarmax = ",Prioralpha$Istarmax,fill=TRUE)
cat("   alphamin = ",Prioralpha$alphamin,fill=TRUE)
cat("   alphamax = ",Prioralpha$alphamax,fill=TRUE)
cat("   power = ",Prioralpha$power,fill=TRUE)
cat(" ",fill=TRUE)
if(drawdelta) 
{
   cat("deltabar",fill=TRUE)
   print(deltabar)
   cat("Ad",fill=TRUE)
   print(Ad)
}
cat(" ",fill=TRUE)
cat("MCMC Parms: ",fill=TRUE)
cat(paste("s=",round(s,3)," w= ",w," R= ",R," keep= ",keep," maxuniq= ",maxuniq,
          " gridsize for lambda hyperparms= ",gridsize),fill=TRUE)
cat("",fill=TRUE)
#
# allocate space for draws
#
if(drawdelta) Deltadraw=matrix(double((floor(R/keep))*nz*nvar),ncol=nz*nvar)
betadraw=array(double((floor(R/keep))*nlgt*nvar),dim=c(nlgt,nvar,floor(R/keep)))
probdraw=matrix(double(floor(R/keep)),ncol=1)
oldbetas=matrix(double(nlgt*nvar),ncol=nvar)
oldll=double(nlgt)
loglike=double(floor(R/keep))
thetaStar=NULL
compdraw=NULL
Istardraw=matrix(double(floor(R/keep)),ncol=1)
alphadraw=matrix(double(floor(R/keep)),ncol=1)
nudraw=matrix(double(floor(R/keep)),ncol=1)
vdraw=matrix(double(floor(R/keep)),ncol=1)
adraw=matrix(double(floor(R/keep)),ncol=1)

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
#
# initialize betas for all units
#
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
# initialize delta
#
if (drawdelta) olddelta=rep(0,nz*nvar)

#
# initialize theta,thetaStar,ind
#
theta=vector("list",nlgt)
for(i in 1:nlgt) {theta[[i]]=list(mu=rep(0,nvar),rooti=diag(nvar))}
ind=double(nlgt)
thetaStar=unique(theta)
nunique=length(thetaStar)
for(j in 1:nunique){
    ind[which(sapply(theta,identical,thetaStar[[j]]))]=j
}
#
# initialize alpha,lambda
#
alpha=1
lambda=list(mubar=rep(0,nvar),Amu=1,nu=nvar+1,V=(nvar+1)*diag(nvar))
#
# fix oldprob (only one comp)
#
oldprob=1

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
      {mgout=rDPGibbs1(oldbetas-Z%*%t(matrix(olddelta,ncol=nz)),theta,thetaStar,ind,
	  lambda,alpha,Prioralpha,lambda_hyper,maxuniq,gridsize)}
   else
      {mgout=rDPGibbs1(oldbetas,theta,thetaStar,ind,
	  lambda,alpha,Prioralpha,lambda_hyper,maxuniq,gridsize)}

   ind=mgout$ind
   lambda=mgout$lambda
   alpha=mgout$alpha
   theta=mgout$theta
   thetaStar=mgout$thetaStar
   Istar=length(thetaStar)
   
   
   # now draw delta | {beta_i}, ind, comps
   if(drawdelta) {olddelta=drawDelta(Z,oldbetas,ind,thetaStar,deltabar,Ad)}
   #
   #  loop over all lgt equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
   #
      for (lgt in 1:nlgt) 
      {
         rootpi=thetaStar[[ind[lgt]]]$rooti
         #  note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
         if(drawdelta) {
            betabar=thetaStar[[ind[lgt]]]$mu+matrix(olddelta,ncol=nz)%*%as.vector(Z[lgt,])}
         else {
            betabar=thetaStar[[ind[lgt]]]$mu }
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
 		  alphadraw[mkeep,]=alpha
        Istardraw[mkeep,]=Istar
        adraw[mkeep,]=lambda$Amu
        nudraw[mkeep,]=lambda$nu
        vdraw[mkeep,]=lambda$V[1,1]/lambda$nu
        loglike[mkeep]=sum(oldll)
        if(drawdelta) Deltadraw[mkeep,]=olddelta
        compdraw[[mkeep]]=list(list(mu=mgout$thetaNp1[[1]],rooti=mgout$thetaNp1[[2]]))
      }
        
}
ctime=proc.time()[3]
cat(" Total Time Elapsed: ",round((ctime-itime)/60,2),fill=TRUE)
if(drawdelta){
   attributes(Deltadraw)$class=c("bayesm.mat","mcmc")
   attributes(Deltadraw)$mcpar=c(1,R,keep)}
attributes(betadraw)$class=c("bayesm.hcoef")
nmix=list(probdraw=probdraw,zdraw=NULL,compdraw=compdraw)
attributes(nmix)$class="bayesm.nmix"
attributes(adraw)$class=c("bayesm.mat","mcmc")
attributes(nudraw)$class=c("bayesm.mat","mcmc")
attributes(vdraw)$class=c("bayesm.mat","mcmc")
attributes(Istardraw)$class=c("bayesm.mat","mcmc")
attributes(alphadraw)$class=c("bayesm.mat","mcmc")
if(drawdelta) 
   {return(list(Deltadraw=Deltadraw,betadraw=betadraw,nmix=nmix,alphadraw=alphadraw,Istardraw=Istardraw,
	            adraw=adraw,nudraw=nudraw,vdraw=vdraw,loglike=loglike))} 
else
   {return(list(betadraw=betadraw,nmix=nmix,alphadraw=alphadraw,Istardraw=Istardraw,
	            adraw=adraw,nudraw=nudraw,vdraw=vdraw,loglike=loglike))} 
}
