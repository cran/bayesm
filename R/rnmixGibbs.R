rnmixGibbs= 
function(Data,Prior,Mcmc)
{
#
# Revision History: 
#   P. Rossi 3/05
#
# purpose: do Gibbs sampling inference for a mixture of multivariate normals
#
# arguments:
#     Data is a list of y which is an n x k matrix of data -- each row
#       is an iid draw from the normal mixture
#     Prior is a list of (Mubar,A,nu,V,a,ncomp)
#       ncomp is required
#       if elements of the prior don't exist, defaults are assumed
#     Mcmc is a list of R and keep (thinning parameter)
# Output:
#     list with elements
#     pdraw -- R/keep x ncomp array of mixture prob draws
#     zdraw -- R/keep x nobs array of indicators of mixture comp identity for each obs
#     compsdraw -- list of R/keep lists of lists of comp parm draws
#        e.g. compsdraw[[i]] is ith draw -- list of ncomp lists
#             compsdraw[[i]][[j]] is list of parms for jth normal component
#             if jcomp=compsdraw[[i]][j]]
#                        ~N(jcomp[[1]],Sigma), Sigma = t(R)%*%R, R^{-1} = jcomp[[2]]
#
# Model:
#        y_i ~ N(mu_ind,Sigma_ind)
#        ind ~ iid multinomial(p)  p is a 1x ncomp vector of probs
# Priors:
#        mu_j ~ N(mubar,Sigma (x) A^-1)
#        mubar=vec(Mubar)
#        Sigma_j ~ IW(nu,V)
#        note: this is the natural conjugate prior -- a special case of multivariate 
#              regression
#        p ~ Dirchlet(a)
#
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
if(missing(Prior)) {pandterm("requires Prior argument <min of Prior$ncomp>")}
else
   {
    if(is.null(Prior$ncomp)) {pandterm("requires number of mix comps -- Prior$ncomp")}
       else {ncomp=Prior$ncomp}
    if(is.null(Prior$Mubar)) {Mubar=matrix(rep(0,dimy),nrow=1)} 
       else {Mubar=Prior$Mubar}
    if(is.null(Prior$A)) {A=matrix(c(.01),ncol=1)} 
       else {A=Prior$A}
    if(is.null(Prior$nu)) {nu=dimy+2} 
       else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(dimy)} 
       else {V=Prior$V}
    if(is.null(Prior$a)) {a=c(rep(5,ncomp))}
       else {a=Prior$a}
   }
#
# check dimensions of Priors
#
if(ncol(A) != nrow(A) || ncol(A) != 1)
   {pandterm(paste("bad dimensions for A",dim(A)))}
if(!is.matrix(Mubar))
   {pandterm("Mubar must be a matrix")}
if(nrow(Mubar) != 1 || ncol(Mubar) != dimy) 
   {pandterm(paste("bad dimensions for Mubar",dim(Mubar)))}
if(ncol(V) != nrow(V) || ncol(V) != dimy)
   {pandterm(paste("bad dimensions for V",dim(V)))}
if(length(a) != ncomp)
   {pandterm(paste("a wrong length, length= ",length(a)))}
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
# print out the problem
#
cat(" Starting Gibbs Sampler for Mixture of Normals",fill=TRUE)
cat(" ",nobs," observations on ",dimy," dimensional data",fill=TRUE)
cat("     using ",ncomp," mixture components",fill=TRUE)
cat(" ",fill=TRUE)
cat(" Prior Parms: ",fill=TRUE)
cat("  mu_j ~ N(mubar,Sigma (x) A^-1)",fill=TRUE)
cat("  mubar = ",fill=TRUE)
print(Mubar)
cat("  precision parm for prior variance of mu vectors (A)= ",A,fill=TRUE)
cat("  Sigma_j ~ IW(nu,V) nu= ",nu,fill=TRUE)
cat("  V =",fill=TRUE)
print(V)
cat("  Dirichlet parameters ",fill=TRUE)
print(a)
cat(" ",fill=TRUE)
cat(" Mcmc Parms: R= ",R," keep= ",keep,fill=TRUE)

pdraw=matrix(double(floor(R/keep)*ncomp),ncol=ncomp)
zdraw=matrix(double(floor(R/keep)*nobs),ncol=nobs)
compdraw=list()
compsd=list()

#
# set initial values of z
#
ninc = floor(nobs/ncomp)
z = c()
for(i in 1:(ncomp-1)) z = c(z,rep(i,ninc))
z = c(z,rep(ncomp,nobs-length(z)))
cat(" ",fill=TRUE)
cat("starting value for z",fill=TRUE)
print(table(z))
cat(" ",fill=TRUE)
p=c(rep(1,ncomp))/ncomp # note this is not used


#
# start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration (est time to end -min) ",fill=TRUE)
fsh()
for(rep in 1:R)
{
   out = rmixGibbs(y,Mubar,A,nu,V,a,p,z,compsd)
   compsd=out$comps
   p=out$p
   z=out$z
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
      pdraw[mkeep,]=p
      zdraw[mkeep,]=z
      compdraw[[rep]]=compsd
      }
}
return(list(probdraw=pdraw,zdraw=zdraw,compdraw=compdraw))
}
