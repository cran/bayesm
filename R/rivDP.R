rivDP = 
function(Data,Prior,Mcmc) 
{
#
# revision history:
#   P. Rossi 1/06
#   added draw of alpha 2/06
#   added automatic scaling 2/06
#   removed reqfun  7/07 -- now functions are in rthetaDP
#   fixed initialization of theta 3/09
#
# purpose: 
#   draw from posterior for linear I.V. model with DP process for errors
#
# Arguments:
#   Data -- list of z,w,x,y
#        y is vector of obs on lhs var in structural equation
#        x is "endogenous" var in structural eqn
#        w is matrix of obs on "exogenous" vars in the structural eqn
#        z is matrix of obs on instruments
#   Prior -- list of md,Ad,mbg,Abg,mubar,Amu,nuV
#        md is prior mean of delta
#        Ad is prior prec
#        mbg is prior mean vector for beta,gamma
#        Abg is prior prec of same
#        lamda is a list of prior parms for DP draw
#              mubar is prior mean of means for "errors"
#              Amu is scale precision parm for means
#              nu,V parms for IW on Sigma (idential priors for each normal comp
#        alpha prior parm for DP process (weight on base measure)
#           or starting value if there is a prior on alpha (requires element Prioralpha)
#        Prioralpha list of hyperparms for draw of alpha (alphamin,alphamax,power,n)
#
#   Mcmc -- list of R,keep,starting values for delta,beta,gamma,theta
#        maxuniq is maximum number of unique theta values
#        R is number of draws
#        keep is thinning parameter
#        SCALE if scale data, def: TRUE
#        gridsize is the gridsize parm for alpha draws
#
#   Output: 
#      list of draws of delta,beta,gamma and thetaNp1 which is used for
#      predictive distribution of errors (density estimation)
# 
#   Model:
#
#    x=z'delta + e1
#    y=beta*x + w'gamma + e2
#        e1,e2 ~ N(theta_i)
#
#   Priors
#   delta ~ N(md,Ad^-1)
#   vec(beta,gamma) ~ N(mbg,Abg^-1)
#   theta ~ DPP(alpha|lambda)
#
#
#   extract data and check dimensios
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of z,w,x,y")}
    if(is.null(Data$w)) isgamma=FALSE else isgamma=TRUE
    if(isgamma) w = Data$w #matrix
    if(is.null(Data$z)) {pandterm("Requires Data element z")}
    z=Data$z
    if(is.null(Data$x)) {pandterm("Requires Data element x")}
    x=as.vector(Data$x)
    if(is.null(Data$y)) {pandterm("Requires Data element y")}
    y=as.vector(Data$y)

#
# check data for validity
#
n=length(y)
if(isgamma)
   {if(!is.matrix(w)) {pandterm("w is not a matrix")}
   dimg=ncol(w)
   if(n != nrow(w) ) {pandterm("length(y) ne nrow(w)")}}

if(!is.matrix(z)) {pandterm("z is not a matrix")}
dimd=ncol(z)
if(n != length(x) ) {pandterm("length(y) ne length(x)")}
if(n != nrow(z) ) {pandterm("length(y) ne nrow(z)")}


#
# extract elements corresponding to the prior
#
if(missing(Prior))
   {
    md=c(rep(0,dimd)) 
    Ad=diag(0.01,dimd) 
    if(isgamma) dimbg=1+dimg else dimbg=1
    mbg=c(rep(0,dimbg)) 
    Abg=diag(0.01,dimbg) 
 

    gamma= .5772156649015328606  
    Istarmin=1
    alphamin=exp(digamma(Istarmin)-log(gamma+log(n)))
    Istarmax=floor(.1*n)
    alphamax=exp(digamma(Istarmax)-log(gamma+log(n)))
    power=.8
    Prioralpha=list(n=n,alphamin=alphamin,alphamax=alphamax,power=power)

    lambda=list(mubar=c(0,0),Amu=.2,nu=3.4,V=1.7*diag(2))
   }

else  
   { 
    if(is.null(Prior$md)) md=c(rep(0,dimd)) else md=Prior$md
    if(is.null(Prior$Ad)) Ad=diag(0.01,dimd) else md=Prior$Ad
    if(isgamma) dimbg=1+dimg else dimbg=1
    if(is.null(Prior$mbg)) mbg=c(rep(0,dimbg)) else md=Prior$mbg
    if(is.null(Prior$Abg)) Abg=diag(0.01,dimbg) else md=Prior$Abg


    if(!is.null(Prior$Prioralpha))
       {Prioralpha=Prior$Prioralpha}
    else
       {gamma= .5772156649015328606  
        Istarmin=1
        alphamin=exp(digamma(Istarmin)-log(gamma+log(n)))
        Istarmax=floor(.1*n)
        alphamax=exp(digamma(Istarmax)-log(gamma+log(n)))
        power=.8
        Prioralpha=list(n=n,alphamin=alphamin,alphamax=alphamax,power=power)}

     if(!is.null(Prior$lambda))
       {lambda=Prior$lambda}
     else
       {lambda=list(mubar=c(0,0),Amu=.2,nu=3.4,V=1.7*diag(2))}
    }

#
# obtain starting values for MCMC
#
# we draw need inital values of delta, theta and indic
#

if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
theta=NULL
if(!is.null(Mcmc$delta)) 
   {delta = Mcmc$delta}
else
   {lmxz = lm(x~z,data.frame(x=x,z=z))
    delta = lmxz$coef[2:(ncol(z)+1)]}
if(!is.null(Mcmc$theta))
  {theta=Mcmc$theta }
else
  {onecomp=list(mu=c(0,0),rooti=diag(2))
   theta=vector("list",length(y))
   for(i in 1:n) {theta[[i]]=onecomp}
   }
dimd = length(delta)
if(is.null(Mcmc$maxuniq))
   {maxuniq=200}
else
   {maxuniq=Mcmc$maxuniq}
if(is.null(Mcmc$R)) {pandterm("requres Mcmc argument, R")}
R = Mcmc$R
if(is.null(Mcmc$keep))
   {keep=1}
else
   {keep=Mcmc$keep}
if(is.null(Mcmc$gridsize))
   {gridsize=20}
else
   {gridsize=Mcmc$gridsize}
if(is.null(Mcmc$SCALE))
  {SCALE=TRUE}
else
  {SCALE=Mcmc$SCALE}


#
# scale and center
#
if(SCALE){
  scaley=sqrt(var(y))
  scalex=sqrt(var(x))
  meany=mean(y)
  meanx=mean(x)
  meanz=apply(z,2,mean)
  y=(y-meany)/scaley; x=(x-meanx)/scalex
  z=scale(z,center=TRUE,scale=FALSE)
  if(isgamma) {meanw=apply(w,2,mean);  w=scale(w,center=TRUE,scale=FALSE)}
}

#
# print out model
#
cat(" ",fill=TRUE)
cat("Starting Gibbs Sampler for Linear IV Model With DP Process Errors",fill=TRUE)
cat(" ",fill=TRUE)
cat(" nobs= ",n,"; ",ncol(z)," instruments",fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("mean of delta ",fill=TRUE)
print(md)
cat(" ",fill=TRUE)
cat("Adelta",fill=TRUE)
print(Ad)
cat(" ",fill=TRUE)
cat("mean of beta/gamma",fill=TRUE)
print(mbg)
cat(" ",fill=TRUE)
cat("Abeta/gamma",fill=TRUE)
print(Abg)
cat(" ",fill=TRUE)
cat("lambda contains: ", fill=TRUE)
cat("mu Prior Parms:",fill=TRUE)
cat("mubar= ",lambda$mubar,fill=TRUE)
cat("Amu= ",lambda$Amu,fill=TRUE)
cat(" ",fill=TRUE)
cat("Sigma Prior Parms:",fill=TRUE)
cat("nu= ",lambda$nu," V=",fill=TRUE)
print(lambda$V)
cat("  ",fill=TRUE)
cat("Parameters of Prior on Dirichlet Process parm (alpha)",fill=TRUE)
cat("alphamin= ",Prioralpha$alphamin," alphamax= ",Prioralpha$alphamax," power=",
        Prioralpha$power,fill=TRUE)
cat("alpha values correspond to Istarmin = ",Istarmin," Istarmax = ",Istarmax,fill=TRUE)
cat(" ",fill=TRUE)
cat("MCMC parms: R= ",R," keep= ",keep,fill=TRUE)
cat("  maximum number of unique thetas= ",maxuniq,fill=TRUE)
cat("  gridsize for alpha draws= ",gridsize,fill=TRUE)
cat("  SCALE data= ",SCALE,fill=TRUE)
cat(" ",fill=TRUE)


#
# define needed functions
#
#
#
# --------------------------------------------------------------------------------------------
#
#
get_ytxt=function(y,z,delta,x,w,ncomp,indic,comps){
yt=NULL; xt=NULL;
if(missing(w)) isw=FALSE else isw=TRUE
if(isw) ncolw=ncol(w)
for (k in 1:ncomp)
{ 
  nobs=sum(indic==k)
  if(nobs > 0) 
     {
     if(isw) wk=matrix(w[indic==k,],ncol=ncolw)
     zk=matrix(z[indic==k,],ncol=length(delta))
     yk=y[indic==k]
     xk=matrix(x[indic==k],ncol=1)
     Sigma=backsolve(comps[[k]][[2]],diag(2))
     Sigma=crossprod(Sigma)
     mu=comps[[k]][[1]]
     e1 = as.vector(xk-zk%*%delta)
     ee2 = mu[2] +(Sigma[1,2]/Sigma[1,1])*(e1-mu[1])
     sig = sqrt(Sigma[2,2]-(Sigma[1,2]^2/Sigma[1,1]))
     yt = c(yt,(yk-ee2)/sig)
     if(isw) 
        {xt = rbind(xt,(cbind(xk,wk)/sig))}
     else
        {xt=rbind(xt,xk/sig)}
     }
}
return(list(xt=xt,yt=yt))
}
#
#
# --------------------------------------------------------------------------------------------
#
#
get_ytxtd=function(y,z,beta,gamma,x,w,ncomp,indic,comps,dimd){
yt=NULL; xtd=NULL;
if(missing(w)) isw=FALSE else isw=TRUE
if(isw) ncolw=ncol(w)
C = matrix(c(1,beta,0,1),nrow=2)
for (k in 1:ncomp)
   {
    nobs=sum(indic==k)
    if(nobs > 0) 
     {
      xtdk=matrix(nrow=2*nobs,ncol=dimd)
      ind=seq(1,(2*nobs-1),by=2)
      if(isw) wk=matrix(w[indic==k,],ncol=ncolw)
      zk=matrix(z[indic==k,],ncol=dimd)
      zveck=as.vector(t(zk))
      yk=y[indic==k]
      xk=x[indic==k]
      Sigma=backsolve(comps[[k]][[2]],diag(2))
      Sigma=crossprod(Sigma)
      mu=comps[[k]][[1]]
      B = C%*%Sigma%*%t(C)
      L = t(chol(B))
      Li=backsolve(L,diag(2),upper.tri=FALSE)
      if(isw) {u=as.vector((yk-wk%*%gamma-mu[2]-beta*mu[1]))}
      else {u=as.vector((yk-mu[2]-beta*mu[1]))}
      ytk = as.vector(Li %*% rbind((xk-mu[1]),u))

      z2=rbind(zveck,beta*zveck)
      z2=Li%*%z2
      zt1=z2[1,]
      zt2=z2[2,]

      dim(zt1)=c(dimd,nobs)
      zt1=t(zt1)
      dim(zt2)=c(dimd,nobs)
      zt2=t(zt2)

      xtdk[ind,]=zt1
      xtdk[-ind,]=zt2

      yt=c(yt,ytk)
      xtd=rbind(xtd,xtdk)
    }
   }
return(list(yt=yt,xtd=xtd))
}
#
#
# --------------------------------------------------------------------------------------------
#
#
rthetaDP= function(maxuniq,alpha,lambda,Prioralpha,theta,thetaStar,indic,q0v,y,gridsize){
# 
#  function to make one draw from DP process 
#
#  P. Rossi 1/06
#  added draw of alpha 2/06
#  removed lambdaD,etaD and function arguments 5/06
#  removed thetaStar argument to .Call and creation of newthetaStar 7/06
#  removed q0 computations as eta is not drawn  7/06
#  changed for new version of thetadraw and removed calculation of thetaStar before
#    .Call  7/07
#
#      y(i) ~ f(y|theta[[i]],eta)
#      theta ~ DP(alpha,G(lambda))
#              note: eta is not used
#output:
#   list with components:
#      thetaDraws: list, [[i]] is a list of the ith draw of the n theta's
#                  where n is the length of the input theta and nrow(y)
#      thetaNp1Draws: list, [[i]] is ith draw of theta_{n+1}
#args:
#   maxuniq: the maximum number of unique thetaStar values -- an error will be raised
#            if this is exceeded
#   alpha,lambda: starting values (or fixed DP prior values if not drawn).
#   Prioralpha: list of hyperparms of alpha prior
#   theta: list of starting value for theta's
#   thetaStar: list of unique values of theta, thetaStar[[i]]
#   indic:  n vector of indicator for which unique theta (in thetaStar)
#   y: is a matrix nxk
#         thetaStar: list of unique values of theta, thetaStar[[i]]
#   q0v:a double vector with the same number of rows as y, giving \Int f(y(i)|theta,eta) dG_{lambda}(theta).
#
#  define needed functions for rthetaDP
# -----------------------------------------------------------------------------------------------
   pandterm = function(message) {
        stop(message, call. = FALSE) }
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
# -----------------------------------------------------------------------------------------------
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
#
# -----------------------------------------------------------------------------------------
#
#
GD=function(lambda){
#
# function to draw from prior for Multivariate Normal Model
#
# mu|Sigma ~ N(mubar,Sigma x Amu^-1)
# Sigma ~ IW(nu,V)
#
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
#
# -------------------------------------------------------------------------------------------
#
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
#  END OF REQUIRED FUNCTIONS AREA
# --------------------------------------------------------------------------------------------
#

   n = length(theta)

   eta=NULL    # note eta is not used
   thetaNp1=NULL

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

   #alpha
   alpha=alphaD(Prioralpha,nunique,gridsize)
   
   return(list(theta=theta,indic=indic,thetaStar=thetaStar,
               thetaNp1=thetaNp1,alpha=alpha,Istar=nunique))
}
#
#
# -----------------------------------------------------------------------------------------
#
#
q0=function(y,lambda,eta){
#
# function to compute a vector of int f(y[i]|theta) p(theta|lambda)dlambda
#     here p(theta|lambda) is G0 the base prior
#
# implemented for a multivariate normal data density and standard conjugate
# prior:
#    theta=list(mu,Sigma)
#    f(y|theta) is N(mu,Sigma)
#    lambda=list(mubar,Amu,nu,V)
#       mu|Sigma ~ N(mubar,Sigma (x) Amu^-1)
#       Sigma ~ IW(nu,V)
#
# arguments:
#    Y is n x k matrix of observations
#    eta is not used
#    lambda=list(mubar,Amu,nu,V)
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
#
#
# --------------------------------------------------------------------------------------------
#
#
#    END OF REQUIRED FUNCTIONS AREA
#
#
#initialize comps,indic,ncomp
comps=unique(theta)
ncomp=length(comps)
indic=double(n)
for(j in 1:ncomp){
      indic[which(sapply(theta,identical,comps[[j]]))]=j
   }
# initialize eta
eta=NULL
#
# initialize alpha
alpha=1

# reserve space for draws
#
deltadraw = matrix(double(floor(R/keep)*dimd),ncol=dimd)
betadraw = rep(0.0,floor(R/keep))
alphadraw=double(floor(R/keep))
Istardraw=double(floor(R/keep))
if(isgamma) gammadraw = matrix(double(floor(R/keep)*dimg),ncol=dimg)
thetaNp1draw=vector("list",R)

#
# start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration (est time to end -min) ",fill=TRUE)
fsh()

for(rep in 1:R) {

   # draw beta and gamma
      if(isgamma) 
         {out=get_ytxt(y=y,z=z,delta=delta,x=x,w=w,
          ncomp=ncomp,indic=indic,comps=comps)}
      else
         {out=get_ytxt(y=y,z=z,delta=delta,x=x,
          ncomp=ncomp,indic=indic,comps=comps)}
         
      bg = breg(out$yt,out$xt,mbg,Abg)  
      beta = bg[1]
      if(isgamma) gamma = bg[2:length(bg)]

   # draw delta
      if(isgamma)
         {out=get_ytxtd(y=y,z=z,beta=beta,gamma=gamma,
          x=x,w=w,ncomp=ncomp,indic=indic,comps=comps,dimd=dimd)}
      else
         {out=get_ytxtd(y=y,z=z,beta=beta,
          x=x,ncomp=ncomp,indic=indic,comps=comps,dimd=dimd)}
	
      delta = breg(out$yt,out$xtd,md,Ad)

    # DP process stuff- theta | lambda
      if(isgamma) {Err = cbind(x-z%*%delta,y-beta*x-w%*%gamma)}
      else {Err = cbind(x-z%*%delta,y-beta*x)}
      q0v = q0(Err,lambda,eta)
      DPout=rthetaDP(maxuniq=maxuniq,alpha=alpha,lambda=lambda,Prioralpha=Prioralpha,theta=theta,
                     thetaStar=comps,indic=indic,q0v=q0v,y=Err,gridsize=gridsize)
      indic=DPout$indic
      theta=DPout$theta
      comps=DPout$thetaStar
      alpha=DPout$alpha
      Istar=DPout$Istar
      ncomp=length(comps)
   
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
      deltadraw[mkeep,]=delta
      betadraw[mkeep]=beta
      alphadraw[mkeep]=alpha
      Istardraw[mkeep]=Istar
      if(isgamma) gammadraw[mkeep,]=gamma
      thetaNp1draw[[mkeep]]=list(DPout$thetaNp1)
      }
}
#
# rescale
#
if(SCALE){
   deltadraw=deltadraw*scalex
   betadraw=betadraw*scaley/scalex
   if(isgamma) {gammadraw=gammadraw*scaley}
}
ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')

nmix=list(probdraw=matrix(c(rep(1,length(thetaNp1draw))),ncol=1),zdraw=NULL,compdraw=thetaNp1draw)
#
# densitymix is in the format to be used with the generic mixture of normals plotting
# methods (plot.bayesm.nmix)
#
attributes(nmix)$class=c("bayesm.nmix")

attributes(deltadraw)$class=c("bayesm.mat","mcmc")
attributes(deltadraw)$mcpar=c(1,R,keep)
attributes(betadraw)$class=c("bayesm.mat","mcmc")
attributes(betadraw)$mcpar=c(1,R,keep)
attributes(alphadraw)$class=c("bayesm.mat","mcmc")
attributes(alphadraw)$mcpar=c(1,R,keep)
attributes(Istardraw)$class=c("bayesm.mat","mcmc")
attributes(Istardraw)$mcpar=c(1,R,keep)
if(isgamma){
   attributes(gammadraw)$class=c("bayesm.mat","mcmc")
   attributes(gammadraw)$mcpar=c(1,R,keep)}

if(isgamma) 
   { return(list(deltadraw=deltadraw,betadraw=betadraw,alphadraw=alphadraw,Istardraw=Istardraw,
                 gammadraw=gammadraw,nmix=nmix))}
   else
   { return(list(deltadraw=deltadraw,betadraw=betadraw,alphadraw=alphadraw,Istardraw=Istardraw,
                 nmix=nmix))}
}


