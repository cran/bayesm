rivGibbs=
function(Data,Prior,Mcmc) 
{
#
# revision history:
#    R. McCulloch original version 2/05 
#    p. rossi 3/05 
#    p. rossi 1/06 -- fixed error in nins
#    p. rossi 1/06 -- fixed def Prior settings for nu,V
#    3/07 added classes
#
# purpose: 
#   draw from posterior for linear I.V. model
#
# Arguments:
#   Data -- list of z,w,x,y
#        y is vector of obs on lhs var in structural equation
#        x is "endogenous" var in structural eqn
#        w is matrix of obs on "exogenous" vars in the structural eqn
#        z is matrix of obs on instruments
#   Prior -- list of md,Ad,mbg,Abg,nu,V
#        md is prior mean of delta
#        Ad is prior prec
#        mbg is prior mean vector for beta,gamma
#        Abg is prior prec of same
#        nu,V parms for IW on Sigma
#
#   Mcmc -- list of R,keep 
#        R is number of draws
#        keep is thinning parameter
#
#   Output: 
#      list of draws of delta,beta,gamma and Sigma
# 
#   Model:
#
#    x=z'delta + e1
#    y=beta*x + w'gamma + e2
#        e1,e2 ~ N(0,Sigma)
#
#   Priors
#   delta ~ N(md,Ad^-1)
#   vec(beta,gamma) ~ N(mbg,Abg^-1)
#   Sigma ~ IW(nu,V)
#
#   check arguments
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of z,w,x,y")}
    if(is.null(Data$z)) {pandterm("Requires Data element z")}
    z=Data$z
    if(is.null(Data$w)) {pandterm("Requires Data element w")}
    w=Data$w
    if(is.null(Data$x)) {pandterm("Requires Data element x")}
    x=Data$x
    if(is.null(Data$y)) {pandterm("Requires Data element y")}
    y=Data$y

#
# check data for validity
#
if(!is.vector(x)) {pandterm("x must be a vector")}
if(!is.vector(y)) {pandterm("y must be a vector")}
n=length(y)
if(!is.matrix(w)) {pandterm("w is not a matrix")}
if(!is.matrix(z)) {pandterm("z is not a matrix")}
dimd=ncol(z)
dimg=ncol(w)
if(n != length(x) ) {pandterm("length(y) ne length(x)")}
if(n != nrow(w) ) {pandterm("length(y) ne nrow(w)")}
if(n != nrow(z) ) {pandterm("length(y) ne nrow(z)")}
#
# check for Prior
#
if(missing(Prior))
   { md=c(rep(0,dimd));Ad=.01*diag(dimd); 
     mbg=c(rep(0,(1+dimg))); Abg=.01*diag((1+dimg));
     nu=3; V=diag(2)}
else
   {
    if(is.null(Prior$md)) {md=c(rep(0,dimd))} 
       else {md=Prior$md}
    if(is.null(Prior$Ad)) {Ad=.01*diag(dimd)} 
       else {Ad=Prior$Ad}
    if(is.null(Prior$mbg)) {mbg=c(rep(0,(1+dimg)))} 
       else {mbg=Prior$mbg}
    if(is.null(Prior$Abg)) {Abg=.01*diag((1+dimg))} 
       else {Abg=Prior$Abg}
    if(is.null(Prior$nu)) {nu=3}
       else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(2)}
       else {V=Prior$V}
   }
#
# check dimensions of Priors
#
if(ncol(Ad) != nrow(Ad) || ncol(Ad) != dimd || nrow(Ad) != dimd) 
   {pandterm(paste("bad dimensions for Ad",dim(Ad)))}
if(length(md) != dimd)
   {pandterm(paste("md wrong length, length= ",length(md)))}
if(ncol(Abg) != nrow(Abg) || ncol(Abg) != (1+dimg) || nrow(Abg) != (1+dimg)) 
   {pandterm(paste("bad dimensions for Abg",dim(Abg)))}
if(length(mbg) != (1+dimg))
   {pandterm(paste("mbg wrong length, length= ",length(mbg)))}
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
# print out model
#
cat(" ",fill=TRUE)
cat("Starting Gibbs Sampler for Linear IV Model",fill=TRUE)
cat(" ",fill=TRUE)
cat(" nobs= ",n,"; ",ncol(z)," instruments; ",ncol(w)," included exog vars",fill=TRUE)
cat("     Note: the numbers above include intercepts if in z or w",fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("mean of delta ",fill=TRUE)
print(md)
cat("Adelta",fill=TRUE)
print(Ad)
cat("mean of beta/gamma",fill=TRUE)
print(mbg)
cat("Abeta/gamma",fill=TRUE)
print(Abg)
cat("Sigma Prior Parms",fill=TRUE)
cat("nu= ",nu," V=",fill=TRUE)
print(V)
cat(" ",fill=TRUE)
cat("MCMC parms: R= ",R," keep= ",keep,fill=TRUE)
cat(" ",fill=TRUE)

deltadraw = matrix(double(floor(R/keep)*dimd),ncol=dimd)
betadraw = rep(0.0,floor(R/keep))
gammadraw = matrix(double(floor(R/keep)*dimg),ncol=dimg)
Sigmadraw = matrix(double(floor(R/keep)*4),ncol=4)

#set initial values
Sigma=diag(2)
delta=c(rep(.1,dimd))

#
# start main iteration loop
#
itime=proc.time()[3]
cat("MCMC Iteration (est time to end -min) ",fill=TRUE)
fsh()
xtd=matrix(nrow=2*n,ncol=dimd)
ind=seq(1,(2*n-1),by=2)
zvec=as.vector(t(z))

for(rep in 1:R) {

    # draw beta,gamma
      e1 = as.vector(x-z%*%delta)
      ee2 = (Sigma[1,2]/Sigma[1,1])*e1
      sig = sqrt(Sigma[2,2]-(Sigma[1,2]^2/Sigma[1,1]))
      yt = (y-ee2)/sig
      xt = cbind(x,w)/sig
      bg = breg(yt,xt,mbg,Abg)  
      beta = bg[1]
      gamma = bg[2:length(bg)]

    # draw delta
      C = matrix(c(1,beta,0,1),nrow=2)
      B = C%*%Sigma%*%t(C)
      L = t(chol(B))
      Li=backsolve(L,diag(2),upper.tri=FALSE)
      u = as.vector((y-w%*%gamma))
      yt = as.vector(Li %*% rbind(x,u))

      z2=rbind(zvec,beta*zvec)
      z2=Li%*%z2
      zt1=z2[1,]
      zt2=z2[2,]
      dim(zt1)=c(dimd,n)
      zt1=t(zt1)
      dim(zt2)=c(dimd,n)
      zt2=t(zt2)
      xtd[ind,]=zt1
      xtd[-ind,]=zt2
      delta = breg(yt,xtd,md,Ad)

    # draw Sigma
      Res = cbind(x-z%*%delta,y-beta*x-w%*%gamma)
      S = crossprod(Res)
      Sigma = rwishart(nu+n,chol2inv(chol(V+S)))$IW
  
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
      gammadraw[mkeep,]=gamma
      Sigmadraw[mkeep,]=Sigma
      }
}

ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')

attributes(deltadraw)$class=c("bayesm.mat","mcmc")
attributes(deltadraw)$mcpar=c(1,R,keep)
attributes(betadraw)$class=c("bayesm.mat","mcmc")
attributes(betadraw)$mcpar=c(1,R,keep)
attributes(gammadraw)$class=c("bayesm.mat","mcmc")
attributes(gammadraw)$mcpar=c(1,R,keep)
attributes(Sigmadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(Sigmadraw)$mcpar=c(1,R,keep)


return(list(deltadraw=deltadraw,betadraw=betadraw,gammadraw=gammadraw,Sigmadraw=Sigmadraw))
}
