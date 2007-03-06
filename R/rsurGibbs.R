rsurGibbs=
function(Data,Prior,Mcmc)
{
# 
# revision history:
#          P. Rossi 9/05
#          3/07 added classes
# Purpose:
#   implement Gibbs Sampler for SUR
# 
# Arguments:
#   Data -- regdata
#           regdata is a list of lists of data for each regression
#           regdata[[i]] contains data for regression equation i
#           regdata[[i]]$y is y, regdata[[i]]$X is X
#           note: each regression can have differing numbers of X vars
#                 but you must have same no of obs in each equation. 
#   Prior -- list of prior hyperparameters
#     betabar,A      prior mean, prior precision
#     nu, V          prior on Sigma
#   Mcmc -- list of MCMC parms
#     R number of draws
#     keep -- thinning parameter
# 
# Output: 
#   list of betadraw,Sigmadraw
#
# Model:
#   y_i = X_ibeta + e_i  
#          y is nobs x 1
#          X is nobs x k_i
#          beta is k_i x 1 vector of coefficients
#          i=1,nreg total regressions
#
#         (e_1,k,...,e_nreg,k) ~ N(0,Sigma) k=1,...,nobs
#
#   we can also write as stacked regression
#   y = Xbeta+e
#       y is nobs*nreg x 1,X is nobs*nreg x (sum(k_i))
#   routine draws beta -- the stacked vector of all coefficients
#
# Priors:  beta ~ N(betabar,A^-1)
#          Sigma ~ IW(nu,V)
# 
#
# check arguments
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata")}
    if(is.null(Data$regdata)) {pandterm("Requires Data element regdata")}
    regdata=Data$regdata
#
# check regdata for validity
#
nreg=length(regdata)
nobs=length(regdata[[1]]$y)
nvar=0
indreg=double(nreg+1)
y=NULL
for (reg in 1:nreg) {
   if(length(regdata[[reg]]$y) != nobs || nrow(regdata[[reg]]$X) != nobs)
      {pandterm(paste("incorrect dimensions for regression",reg))}
   else
      {indreg[reg]=nvar+1
       nvar=nvar+ncol(regdata[[reg]]$X); y=c(y,regdata[[reg]]$y)}
} 
indreg[nreg+1]=nvar+1
#
# check for Prior
#
if(missing(Prior))
   { betabar=c(rep(0,nvar)); A=.01*diag(nvar); nu=nreg+3; V=nu*diag(nreg)}
else
   {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))} 
       else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=.01*diag(nvar)} 
       else {A=Prior$A}
    if(is.null(Prior$nu)) {nu=nreg+3}
       else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(nreg)}
       else {ssq=Prior$V}
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
    if(is.null(Mcmc$R)) {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
   }
#
# print out problem
#
cat(" ", fill=TRUE)
cat("Starting Gibbs Sampler for SUR Regression Model",fill=TRUE)
cat("  with ",nreg," regressions",fill=TRUE)
cat("  and  ",nobs," observations for each regression",fill=TRUE)
cat(" ", fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat("nu = ",nu,fill=TRUE)
cat("V = ",fill=TRUE)
print(V)
cat(" ", fill=TRUE)
cat("MCMC parms: ",fill=TRUE)
cat("R= ",R," keep= ",keep,fill=TRUE)
cat(" ",fill=TRUE)

Sigmadraw=matrix(double(floor(R*nreg*nreg/keep)),ncol=nreg*nreg)
betadraw=matrix(double(floor(R*nvar/keep)),ncol=nvar)


#
# set initial value of Sigma
#
E=matrix(double(nobs*nreg),ncol=nreg)
for (reg in 1:nreg) {
    E[,reg]=lm(y~.-1,data=data.frame(y=regdata[[reg]]$y,regdata[[reg]]$X))$residuals
}
Sigma=crossprod(E)/nobs
L=t(backsolve(chol(Sigma),diag(nreg)))
Y=y
dim(Y)=c(nobs,nreg)
Xti=matrix(0,ncol=nvar,nrow=nreg*nobs)

itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
fsh()

for (rep in 1:R)
{
#
#     first draw beta | Sigma
#
# compute Xtilde
#
  for (reg in 1:nreg){
     Xti[,indreg[reg]:(indreg[reg+1]-1)]=L[,reg]%x%regdata[[reg]]$X
  }
  IR=backsolve(chol(crossprod(Xti)+A),diag(nvar))
#
# compute ytilde
  yti=as.vector(Y%*%t(L))
  btilde=crossprod(t(IR))%*%(crossprod(Xti,yti)+A%*%betabar)
  beta = btilde + IR%*%rnorm(nvar)
#
#    now draw Sigma | beta
#
  for(reg in 1:nreg){
     E[,reg]=regdata[[reg]]$y-regdata[[reg]]$X%*%beta[indreg[reg]:(indreg[reg+1]-1)]
  }
  Sigma=rwishart(nu+nobs,chol2inv(chol(crossprod(E)+V)))$IW
  L=t(backsolve(chol(Sigma),diag(nreg)))
#
#       print time to completion and draw # every 100th draw
#
  if(rep%%100 == 0)
    {ctime=proc.time()[3]
    timetoend=((ctime-itime)/rep)*(R-rep)
    cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
    fsh()}

  if(rep%%keep == 0) 
    {mkeep=rep/keep; betadraw[mkeep,]=beta; Sigmadraw[mkeep,]=Sigma}
}
ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')

attributes(betadraw)$class=c("bayesm.mat","mcmc")
attributes(betadraw)$mcpar=c(1,R,keep)
attributes(Sigmadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(Sigmadraw)$mcpar=c(1,R,keep)

return(list(betadraw=betadraw,Sigmadraw=Sigmadraw))
}
