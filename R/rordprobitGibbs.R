rordprobitGibbs=
 function(Data,Prior,Mcmc)
{
#
# revision history:
#   3/07  Hsiu-Wen Liu
#    
# purpose: 
#   draw from posterior for ordered probit using Gibbs Sampler
#   and metropolis RW
#
# Arguments:
#   Data - list of X,y,k  
#     X is nobs x nvar, y is nobs vector of 1,2,.,k (ordinal variable)
#   Prior - list of A, betabar
#     A is nvar x nvar prior preci matrix
#     betabar is nvar x 1 prior mean
#     Ad is ndstar x ndstar prior preci matrix of dstar (ncut is number of cut-offs being estimated)
#     dstarbar is ndstar x 1 prior mean of dstar
#   Mcmc
#     R is number of draws
#     keep is thinning parameter
#     s is scale parameter of random work Metropolis    
#      
# Output:
#   list of betadraws and cutdraws
#
# Model: 
#    z=Xbeta + e  < 0  e ~N(0,1)
#    y=1,..,k, if z~c(c[k], c[k+1])
#
#    cutoffs = c[1],..,c[k+1]
#    dstar = dstar[1],dstar[k-2]
#    set c[1]=-100, c[2]=0, ...,c[k+1]=100
#
#    c[3]=exp(dstar[1]),c[4]=c[3]+exp(dstar[2]),...,
#    c[k]=c[k-1]+exp(datsr[k-2])
#    
# Note: 1. length of dstar = length of cutoffs - 3
#       2. Be careful in assessing prior parameter, Ad.  .1 is too small for many applications.
#
# Prior: beta ~ N(betabar,A^-1)
#        dstar ~ N(dstarbar, Ad^-1)
#
#
# ----------------------------------------------------------------------
# define functions needed
#
breg1=
function(root,X,y,Abetabar) 
{
#
# p.rossi 12/04
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
# Prior:  beta ~ N(betabar,A^-1)
#
cov=crossprod(root,root)
betatilde=cov%*%(crossprod(X,y)+Abetabar)
betatilde+t(root)%*%rnorm(length(betatilde))
}

#  
#  dstartoc is a fuction to transfer dstar to its cut-off value    

    dstartoc=function(dstar) {c(-100, 0, cumsum(exp(dstar)), 100)} 

# compute conditional likelihood of data given cut-offs
#
   lldstar=function(dstar,y,mu){
           gamma=dstartoc(dstar)
           arg = pnorm(gamma[y+1]-mu)-pnorm(gamma[y]-mu)
           epsilon=1.0e-50
           arg=ifelse(arg < epsilon,epsilon,arg)
           return(sum(log(arg)))
           }


dstarRwMetrop=
function(y,mu,olddstar,s,inc.root,dstarbar,oldll,rootdi){ 
#
# function to execute rw metropolis for the dstar
# y is n vector with element = 1,...,j 
# X is n x k matrix of x values 
# RW increments are N(0,s^2*t(inc.root)%*%inc.root)
# prior on dstar is N(dstarbar,Sigma)  Sigma^-1=rootdi*t(rootdi)
#	inc.root, rootdi are upper triangular
#	this means that we are using the UL decomp of Sigma^-1 for prior 
# olddstar is the current
     
     stay=0   
     dstarc=olddstar + s*t(inc.root)%*%(matrix(rnorm(ncut),ncol=1))
     cll=lldstar(dstarc,y,mu)
     clpost=cll+lndMvn(dstarc,dstarbar,rootdi)
     ldiff=clpost-oldll-lndMvn(olddstar,dstarbar,rootdi)
     alpha=min(1,exp(ldiff))

     if(alpha < 1) {unif=runif(1)} else {unif=0}
     if (unif <= alpha)
          {dstardraw=dstarc; oldll=cll}
     else 
          {dstardraw=olddstar; stay=1}

return(list(dstardraw=dstardraw,oldll=oldll, stay=stay))
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
    if(is.null(Data$k)) {pandterm("Requires Data element k")}
    k=Data$k

nvar=ncol(X)
nobs=length(y)  
ndstar = k-2         # number of dstar being estimated
ncuts = k+1          # number of cut-offs (including zero and two ends)
ncut = ncuts-3       # number of cut-offs being estimated c[1]=-100, c[2]=0, c[k+1]=100 

#
# check data for validity
#
if(length(y) != nrow(X) ) {pandterm("y and X not of same row dim")}
if(  sum(unique(y) %in% (1:k) ) < length(unique(y)) )
  {pandterm("some value of y is not vaild")}

#

#
# check for Prior
#
if(missing(Prior))
   { betabar=c(rep(0,nvar)); A=.01*diag(nvar); Ad=diag(ndstar); dstarbar=c(rep(0,ndstar))}
else
   {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))} 
       else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=.01*diag(nvar)} 
       else {A=Prior$A}
    if(is.null(Prior$Ad)) {Ad=diag(ndstar)} 
       else {Ad=Prior$Ad}
    if(is.null(Prior$dstarbar)) {dstarbar=c(rep(0,ndstar))} 
       else {dstarbar=Prior$dstarbar}
   }
#
# check dimensions of Priors
#

if(ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) 
   {pandterm(paste("bad dimensions for A",dim(A)))}
if(length(betabar) != nvar)
   {pandterm(paste("betabar wrong length, length= ",length(betabar)))}
if(ncol(Ad) != nrow(Ad) || ncol(Ad) != ndstar || nrow(Ad) != ndstar) 
   {pandterm(paste("bad dimensions for Ad",dim(Ad)))}
if(length(dstarbar) != ndstar)
   {pandterm(paste("dstarbar wrong length, length= ",length(dstarbar)))}

#
# check MCMC argument
#
if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
else
   {
    if(is.null(Mcmc$R)) 
       {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
    if(is.null(Mcmc$s)) {s=2.93/sqrt(ndstar)} else {s=Mcmc$s} 
    }
#
# print out problem
#
cat(" ", fill=TRUE)
cat("Starting Gibbs Sampler for Ordered Probit Model",fill=TRUE)
cat("   with ",nobs,"observations",fill=TRUE)
cat(" ", fill=TRUE)
cat("Table of y values",fill=TRUE)
print(table(y))
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat(" ", fill=TRUE)
cat("A",fill=TRUE)
print(A)
cat(" ", fill=TRUE)
cat("dstarbar",fill=TRUE)
print(dstarbar)
cat(" ", fill=TRUE)
cat("Ad",fill=TRUE)
print(Ad)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep,"s= ",s, fill=TRUE) 
cat(" ",fill=TRUE)

betadraw=matrix(double(floor(R/keep)*nvar),ncol=nvar)
cutdraw=matrix(double(floor(R/keep)*ncuts),ncol=ncuts)
dstardraws=matrix(double(floor(R/keep)*ndstar),ncol=ndstar)
staydraw=array(0,dim=c(R/keep))
dstardraw=c(rep(0,ndstar))

sigma=c(rep(1,nrow(X)))
root=chol(chol2inv(chol((crossprod(X,X)+A))))
Abetabar=crossprod(A,betabar)
rootdi=chol(chol2inv(chol(Ad)))

# use (-Hessian+Ad)^(-1) evaluated at betahat as the basis of the 
# covariance matrix for the random walk Metropolis increments 
    
    betahat = chol2inv(chol(crossprod(X,X)))%*% crossprod(X,y)
    dstarini = c(cumsum(c( rep(0.1, ndstar))))     # set initial value for dstar   
    dstarout = optim(dstarini, lldstar, method = "BFGS", hessian=T,
                control = list(fnscale = -1,maxit=500,
                reltol = 1e-06, trace=0), mu=X%*%betahat, y=y)             
    inc.root=chol(chol2inv(chol((-dstarout$hessian+Ad))))  # chol((H+Ad)^-1) 

# set initial values for MCMC  
    
    olddstar = c(rep(0,ndstar))
    beta = betahat    
    cutoffs = dstartoc (olddstar)  
    oldll = lldstar(olddstar,y,mu=X%*%betahat)
 
#
#
#	start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
fsh()
#    print time to completion and draw # every 100th draw
#
for (rep in 1:R) 
{
   # draw z given beta(i-1), sigma, y, cut-offs     
      z = rtrun (X%*%beta, sigma=sigma, a=cutoffs[y] , b=cutoffs[y+1])
  
   # draw beta given z and rest
      beta= breg1(root,X,z, Abetabar)     
   
    # draw gamma given z
      metropout = dstarRwMetrop(y,X%*%beta,olddstar,s,inc.root,dstarbar,oldll,rootdi)     
      olddstar = metropout$dstardraw
      oldll =  metropout$oldll
      cutoffs = dstartoc (olddstar) 
      stay = metropout$stay  


#    print time to completion and draw # every 100th draw

  if(rep%%1000 == 0)
    {ctime=proc.time()[3]
    timetoend=((ctime-itime)/rep)*(R-rep)
    cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
    fsh()}

  if(rep%%keep == 0) 
    {mkeep=rep/keep; cutdraw[mkeep,]=cutoffs; dstardraws[mkeep,]=olddstar;betadraw[mkeep,]=beta;staydraw[mkeep]=stay }
                
}
    accept=1-sum(staydraw)/(R/keep)

ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')

cutdraw=cutdraw[,2:k]
attributes(cutdraw)$class="bayesm.mat"
attributes(betadraw)$class="bayesm.mat"
attributes(dstardraw)$class="bayesm.mat"
attributes(cutdraw)$mcpar=c(1,R,keep)
attributes(betadraw)$mcpar=c(1,R,keep)
attributes(betadraw)$mcpar=c(1,R,keep)

return(list(cutdraw=cutdraw,betadraw=betadraw, dstardraws=dstardraws, accept=accept))
}
