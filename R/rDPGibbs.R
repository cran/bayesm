rDPGibbs= 
function(Prior,Data,Mcmc)
{
#
# Revision History: 
#   5/06 add rthetaDP
#   7/06 include rthetaDP in main body to avoid copy overhead
#   1/08 add scaling
#   2/08 add draw of lambda
#   3/08 changed nu prior support to dim(y) + exp(unif gird on nulim[1],nulim[2])
#
# purpose: do Gibbs sampling for density estimation using Dirichlet process model
#
# arguments:
#     Data is a list of y which is an n x k matrix of data
#     Prior is a list of (alpha,lambda,Prioralpha)
#       alpha: starting value
#       lambda_hyper: hyperparms of prior on lambda
#       Prioralpha: hyperparms of alpha prior; a list of (Istarmin,Istarmax,power)
#       if elements of the prior don't exist, defaults are assumed
#     Mcmc is a list of (R,keep,maxuniq)
#       R: number of draws
#       keep: thinning parameter
#       maxuniq: the maximum number of unique thetaStar values
#
# Output:
#     list with elements
#     alphadraw: vector of length R/keep, [i] is ith draw of alpha
#     Istardraw: vector of length R/keep, [i] is the number of unique theta's drawn from ith iteration
#     adraw
#     nudraw
#     vdraw
#     thetaNp1draws: list, [[i]] is ith draw of theta_{n+1}
#     inddraw: R x n matrix, [,i] is indicators of identity for each obs in ith iteration
#
# Model:
#        y_i ~ f(y|thetai)
#        thetai|G ~ G
#        G|lambda,alpha ~ DP(G|G0(lambda),alpha)
#
# Priors:
#        alpha: starting value
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
#
#
# define needed functions
#
# -----------------------------------------------------------------------------------------
#
q0=function(y,lambda,eta){
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
   rmultinomF=
      function(p) {
       return(sum(runif(1) > cumsum(p))+1)
   }
# -----------------------------------------------------------------------------------------------
   alphaD=function(Prioralpha,Istar,gridsize){
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
yden=function(thetaStar,y,eta){
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
GD=function(lambda){
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
thetaD=function(y,lambda,eta){
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
lndMvn=function (x, mu, rooti) 
{
    return(-(length(x)/2) * log(2 * pi) - 0.5 * sum(((x-mu)%*%rooti)**2) + sum(log(diag(rooti))))
}
# -----------------------------------------------------------------------------------------
   lambdaD=function(lambda,thetastar,alim=c(.01,2),nulim=c(.01,2),vlim=c(.1,5),gridsize=20){
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
   rmultinomF=
      function(p) {
       return(sum(runif(1) > cumsum(p))+1)
   }
echo=function(lst){return(t(lst[[2]]))}
rootiz=function(lst){crossprod(lst[[2]],lst[[1]])}
#
# ------------------------------------------------------------------------------------------

   d=length(thetastar[[1]]$mu)
   Istar=length(thetastar)
   aseq=seq(from=alim[1],to=alim[2],len=gridsize)
   nuseq=d-1+exp(seq(from=nulim[1],to=nulim[2],len=gridsize)) # log uniform grid
   vseq=seq(from=vlim[1],to=vlim[2],len=gridsize)
#
#    "brute" force approach would simply loop over the 
#         "observations" (theta_j) and use log of the appropriate densities.  To vectorize, we
#         notice that the "data" comes via various statistics:
#         1. sum of log(diag(rooti_j)
#         2. sum of tr(V%*%rooti_j%*%t(rooti_j)) where V=vI_d
#         3. quadratic form t(mu_j-0)%*%rooti%*%t(rooti)%*%(mu_j-0)
#     thus, we will compute these first.
#     for documentation purposes, we leave brute force code in comment fields
#
# extract needed info from thetastar list
#
   out=double(Istar*d*d)
   out=sapply(thetastar,echo)
   dim(out)=c(d,Istar*d) # out has the rootis in form: [t(rooti_1), t(rooti_2), ...,t(rooti_Istar)]
   sumdiagriri=sum(colSums(out^2)) #  sum_j tr(rooti_j%*%t(rooti_j))
#   now get diagonals of rooti
   ind=cbind(c(1:(d*Istar)),rep((1:d),Istar))
   out=t(out)
   sumlogdiag=sum(log(out[ind]))
   rimu=sapply(thetastar,rootiz) # columns of rimu contain t(rooti_j)%*%mu_j
   dim(rimu)=c(d,Istar)
   sumquads=sum(colSums(rimu^2)) 
#  
#  draw a  (conditionally indep of nu,v given theta_j)
   lnprob=double(length(aseq))
       #for(i in seq(along=aseq)){
       #for(j in seq(along=thetastar)){
       #lnprob[i]=lnprob[i]+lndMvn(thetastar[[j]]$mu,c(rep(0,d)),thetastar[[j]]$rooti*sqrt(aseq[i]))}
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
       #for(i in seq(along=nuseq)){
       #for(j in seq(along=thetastar)){
       #Sigma_j=crossprod(backsolve(thetastar[[j]]$rooti,diag(d)))
       #lnprob[i]=lnprob[i]+lndIWishart(nuseq[i],V,Sigma_j)}
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
       #for(i in seq(along=vseq)){
       #V=vseq[i]*diag(d)
       #for(j in seq(along=thetastar)){
       #Sigma_j=crossprod(backsolve(thetastar[[j]]$rooti,diag(d)))
       #lnprob[i]=lnprob[i]+lndIWishart(nudraw,V,Sigma_j)}
#    lnprob=Istar*nudraw*d*log(sqrt(vseq))-.5*sumdiagriri*vseq
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
# -----------------------------------------------------------------------------------------

#  check arguments
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of y")}
    if(is.null(Data$y)) {pandterm("Requires Data element y")}
    y=Data$y
#
# check data for validity
#
if(!is.matrix(y)) {pandterm("y must be a matrix")}
nobs=nrow(y)
dimy=ncol(y)
#
# check for Prior
#
alimdef=c(.01,10)
nulimdef=c(.01,3)
vlimdef=c(.1,4)
if(missing(Prior)) {pandterm("requires Prior argument ")}
else
   {
    if(is.null(Prior$lambda_hyper)) {lambda_hyper=list(alim=alimdef,nulim=nulimdef,vlim=vlimdef)}
    else {lambda_hyper=Prior$lambda_hyper;
       if(is.null(lambda_hyper$alim)) {lambda_hyper$alim=alimdef}
       if(is.null(lambda_hyper$nulim)) {lambda_hyper$nulim=nulimdef} 
       if(is.null(lambda_hyper$vlim)) {lambda_hyper$vlim=vlimdef}
       }
    if(is.null(Prior$Prioralpha)) {Prioralpha=list(Istarmin=1,Istarmax=min(50,0.1*nobs),power=0.8)}
    else {Prioralpha=Prior$Prioralpha;
       if(is.null(Prioralpha$Istarmin)) {Prioralpha$Istarmin=1} else {Prioralpha$Istarmin=Prioralpha$Istarmin}
       if(is.null(Prioralpha$Istarmax)) 
             {Prioralpha$Istarmax=min(50,0.1*nobs)} else {Prioralpha$Istarmax=Prioralpha$Istarmax}
       if(is.null(Prioralpha$power)) {Prioralpha$power=0.8}
       }
   }
gamma= .5772156649015328606
Prioralpha$alphamin=exp(digamma(Prioralpha$Istarmin)-log(gamma+log(nobs)))
Prioralpha$alphamax=exp(digamma(Prioralpha$Istarmax)-log(gamma+log(nobs)))
Prioralpha$n=nobs
#
# check Prior arguments for valdity
#
if(lambda_hyper$alim[1]<0) {pandterm("alim[1] must be >0")}
if(lambda_hyper$nulim[1]<0) {pandterm("nulim[1] must be >0")}
if(lambda_hyper$vlim[1]<0) {pandterm("vlim[1] must be >0")}
if(Prioralpha$Istarmin <1){pandterm("Prioralpha$Istarmin must be >= 1")}
if(Prioralpha$Istarmax <= Prioralpha$Istarmin){pandterm("Prioralpha$Istarmin must be > Prioralpha$Istarmax")}
#
# check MCMC argument
#
if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
else
   {
    if(is.null(Mcmc$R)) 
       {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
    if(is.null(Mcmc$maxuniq)) {maxuniq=200} else {maxuniq=Mcmc$maxuniq}
    if(is.null(Mcmc$SCALE)) {SCALE=TRUE} else {SCALE=Mcmc$SCALE}
    if(is.null(Mcmc$gridsize)) {gridsize=20} else {gridsize=Mcmc$gridsize}
   }

#
# print out the problem
#
cat(" Starting Gibbs Sampler for Density Estimation Using Dirichlet Process Model",fill=TRUE)
cat(" ",nobs," observations on ",dimy," dimensional data",fill=TRUE)
cat(" ",fill=TRUE)
cat(" SCALE=",SCALE,fill=TRUE)
cat(" ",fill=TRUE)
cat(" Prior Parms: ",fill=TRUE)
cat("  G0 ~ N(mubar,Sigma (x) Amu^-1)",fill=TRUE)
cat("   mubar = ",0,fill=TRUE)
cat("   Sigma ~ IW(nu,nu*v*I)",fill=TRUE)
cat("   Amu ~ uniform[",lambda_hyper$alim[1],",",lambda_hyper$alim[2],"]",fill=TRUE)
cat("   nu ~ uniform on log grid on [",dimy-1+exp(lambda_hyper$nulim[1]),
             ",",dimy-1+exp(lambda_hyper$nulim[2]),"]",fill=TRUE)
cat("   v ~ uniform[",lambda_hyper$vlim[1],",",lambda_hyper$vlim[2],"]",fill=TRUE)
cat(" ",fill=TRUE)
cat("  alpha ~ (1-(alpha-alphamin)/(alphamax-alphamin))^power",fill=TRUE)
cat("   Istarmin = ",Prioralpha$Istarmin,fill=TRUE)
cat("   Istarmax = ",Prioralpha$Istarmax,fill=TRUE)
cat("   alphamin = ",Prioralpha$alphamin,fill=TRUE)
cat("   alphamax = ",Prioralpha$alphamax,fill=TRUE)
cat("   power = ",Prioralpha$power,fill=TRUE)
cat(" ",fill=TRUE)
cat(" Mcmc Parms: R= ",R," keep= ",keep," maxuniq= ",maxuniq," gridsize for lambda hyperparms= ",gridsize,
        fill=TRUE)
cat(" ",fill=TRUE)

# initialize theta, thetastar, indic

theta=vector("list",nobs)
for(i in 1:nobs) {theta[[i]]=list(mu=rep(0,dimy),rooti=diag(dimy))}
indic=double(nobs)
thetaStar=unique(theta)
nunique=length(thetaStar)
for(j in 1:nunique){
    indic[which(sapply(theta,identical,thetaStar[[j]]))]=j
}
#
# initialize lambda
#
lambda=list(mubar=rep(0,dimy),Amu=.1,nu=dimy+1,V=(dimy+1)*diag(dimy))

#
# initialize alpha
#
alpha=1

alphadraw=double(floor(R/keep))
Istardraw=double(floor(R/keep))
adraw=double(floor(R/keep))
nudraw=double(floor(R/keep))
vdraw=double(floor(R/keep))
thetaNp1draw=vector("list",floor(R/keep))
inddraw=matrix(double((floor(R/keep))*nobs),ncol=nobs)

#
# do scaling
#
if(SCALE){
  dvec=sqrt(apply(y,2,var))
  ybar=apply(y,2,mean)
  y=scale(y,center=ybar,scale=dvec)
  dvec=1/dvec  # R function scale divides by scale
} 
#
# note on scaling
#
#   we model scaled y, z_i=D(y_i-ybar)   D=diag(1/sigma1, ..., 1/sigma_dimy)
#
#   if p_z= 1/R sum(phi(z|mu,Sigma))
#      p_y=1/R sum(phi(y|D^-1mu+ybar,D^-1SigmaD^-1)
#      rooti_y=Drooti_z
#
#   you might want to use quantiles instead, like median and (10,90)
#

#
# start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration (est time to end -min) ",fill=TRUE)
fsh()

for(rep in 1:R)
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

   if(rep%%100==0)
     {
      ctime=proc.time()[3]
      timetoend=((ctime-itime)/rep)*(R-rep)
      cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
      fsh()
      }
   if(rep%%keep ==0)
     {
      mkeep=rep/keep
      alphadraw[mkeep]=alpha
      Istardraw[mkeep]=nunique
      adraw[mkeep]=lambda$Amu
      nudraw[mkeep]=lambda$nu
      vdraw[mkeep]=lambda$V[1,1]/lambda$nu
      if(SCALE){
        thetaNp1[[1]]=thetaNp1[[1]]/dvec+ybar
        if(ncol(y)>1) 
          {thetaNp1[[2]]=diag(dvec)%*%thetaNp1[[2]]}
        else
          {thetaNp1[[2]]=dvec*thetaNp1[[2]]}
      }
      thetaNp1draw[[mkeep]]=list(list(mu=thetaNp1[[1]],rooti=thetaNp1[[2]]))
                            #  here we put the draws into the list of lists of list format useful for
                            #  finite mixture of normals utilities
      inddraw[mkeep,]=indic
      }
}
ctime=proc.time()[3]
cat("Total Time Elapsed= ",round((ctime-itime)/60,2),fill=TRUE)
nmix=list(probdraw=matrix(c(rep(1,nrow(inddraw))),ncol=1),zdraw=inddraw,compdraw=thetaNp1draw)
attributes(nmix)$class="bayesm.nmix"
attributes(alphadraw)$class=c("bayesm.mat","mcmc")
attributes(Istardraw)$class=c("bayesm.mat","mcmc")
attributes(adraw)$class=c("bayesm.mat","mcmc")
attributes(nudraw)$class=c("bayesm.mat","mcmc")
attributes(vdraw)$class=c("bayesm.mat","mcmc")
return(list(alphadraw=alphadraw,Istardraw=Istardraw,adraw=adraw,nudraw=nudraw,
         vdraw=vdraw,nmix=nmix))
}
