rmnpGibbs=
function(Data,Prior,Mcmc) 
{
#
# Revision History:
#   modified by rossi 12/18/04 to include error checking
#   3/07 added classes
#
# purpose:  Gibbs MNP model with full covariance matrix
#
# Arguments:
#   Data contains 
#      p the number of choice alternatives
#      y -- a vector of length n with choices (takes on values from 1, .., p)
#      X -- n(p-1) x k matrix of covariates (including intercepts)
#           note: X is the differenced matrix unlike MNL X=stack(X_1,..,X_n) 
#                 each X_i is (p-1) x nvar
#
#   Prior contains a list of (betabar, A, nu, V)
#      if elements of prior do not exist, defaults are used
#
#   Mcmc is a list of (beta0,sigma0,R,keep)  
#     beta0,sigma0 are intial values, if not supplied defaults are used
#     R is number of draws
#     keep is thinning parm, keep every keepth draw
#
# Output: a list of every keepth betadraw and sigmsdraw
#
#  model: 
#    w_i = X_ibeta + e    e~N(0,Sigma)     note w_i,e are (p-1) x 1
#    y_i = j  if w_ij > w_i-j  j=1,...,p-1
#    y_i = p  if all w_i < 0
#  
#  priors:
#    beta ~ N(betabar,A^-1)
#    Sigma ~ IW(nu,V)
#
#  Check arguments
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of p, y, X")}
  if(is.null(Data$p)) {pandterm("Requires Data element p -- number of alternatives")}
  p=Data$p
  if(is.null(Data$y)) {pandterm("Requires Data element y -- number of alternatives")}
  y=Data$y
  if(is.null(Data$X)) {pandterm("Requires Data element X -- matrix of covariates")}
  X=Data$X
#
# check data for validity
#
levely=as.numeric(levels(as.factor(y)))
if(length(levely) != p) {pandterm(paste("y takes on ",length(levely),
  " values -- must be ",p))}
  bady=FALSE
  for (i in 1:p) 
  {
      if(levely[i] != i) bady=TRUE
  }
cat("Table of y values",fill=TRUE)
print(table(y))
if (bady) {pandterm("Invalid y")}
n=length(y)
k=ncol(X)
pm1=p-1
if(nrow(X)/n != pm1) {pandterm(paste("X has ",nrow(X)," rows; must be = (p-1)n"))}
#
# check for prior elements
#
if(missing(Prior)) 
  { betabar=rep(0,k) ; A=.01*diag(k) ; nu=pm1+3; V=nu*diag(pm1)}
else 
  {if(is.null(Prior$betabar)) {betabar=rep(0,k)} else {betabar=Prior$betabar}
   if(is.null(Prior$A)) {A=.01*diag(k)} else {A=Prior$A}
   if(is.null(Prior$nu)) {nu=pm1+3} else {nu=Prior$nu}
   if(is.null(Prior$V)) {V=nu*diag(pm1)} else {V=Prior$V}}
if(length(betabar) != k) pandterm("length betabar ne k")
if(sum(dim(A)==c(k,k)) != 2) pandterm("A is of incorrect dimension")
if(nu < 1) pandterm("invalid nu value")
if(sum(dim(V)==c(pm1,pm1)) != 2) pandterm("V is of incorrect dimension")
#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R must be included")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$beta0)) {beta0=rep(0,k)} else {beta0=Mcmc$beta0}
if(is.null(Mcmc$sigma0)) {sigma0=diag(pm1)} else {sigma0=Mcmc$sigma0}
if(length(beta0) != k) pandterm("beta0 is not of length k")
if(sum(dim(sigma0) == c(pm1,pm1)) != 2) pandterm("sigma0 is of incorrect dimension")
if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting Gibbs Sampler for MNP",fill=TRUE)
cat("  ",n," obs; ",p," choice alternatives; ",k," indep vars (including intercepts)",fill=TRUE)
cat("  ",R," reps; keeping every ",keep,"th draw",fill=TRUE)
cat(" ",fill=TRUE)
cat("Table of y values",fill=TRUE)
print(table(y))
cat("Prior Parms:",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat("nu",fill=TRUE)
print(nu)
cat("V",fill=TRUE)
print(V)
cat(" ",fill=TRUE)
cat("MCMC Parms:",fill=TRUE)
cat("R= ",R,fill=TRUE)
cat("initial beta= ",beta0,fill=TRUE)
cat("initial sigma= ",sigma0,fill=TRUE)
cat(" ",fill=TRUE)
#
# allocate space for draws
#
sigmadraw=matrix(double(floor(R/keep)*pm1*pm1),ncol=pm1*pm1)
betadraw=matrix(double(floor(R/keep)*k),ncol=k)
wnew=double(nrow(X))
betanew=double(k)

#
#  set initial values of w,beta, sigma (or root of inv)
#
wold=c(rep(0,nrow(X)))
betaold=beta0
C=chol(solve(sigma0))
#
#  C is upper triangular root of sigma^-1 (G) = C'C
#
#  create functions needed
#
drawwc=function(w,mu,y,sigi) {
      .C("draww",w=as.double(w),as.double(mu),as.double(sigi),
        as.integer(length(y)),as.integer(ncol(sigi)),as.integer(y))$w}

draww=
function(w,X,y,beta,sigmai){
#
#   draw latent vector
#
#  	w is n x (p-1) vector
#       X ix n(p-1) x k  matrix
#       y is multinomial 1,..., p
#       beta is k x 1 vector
#       sigmai is (p-1) x (p-1) 
#

Xbeta=as.vector(X%*%beta)
drawwc(w,Xbeta,y,sigmai)
}

itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
for (rep in 1:R) 
   {
   #
   # draw w given beta(rep-1),sigma(rep-1)
   #
   sigmai=crossprod(C)
   wnew=draww(wold,X,y,betaold,sigmai)
   #
   # draw beta given w(rep) and sigma(rep-1)
   #
   #  note:  if Sigma^-1 (G) = C'C then Var(Ce)=CSigmaC' = I
   #  first, transform w_i = X_ibeta + e_i by premultiply by C
   #
   zmat=matrix(cbind(wnew,X),nrow=pm1)
   zmat=C%*%zmat
   zmat=matrix(zmat,nrow=nrow(X))
   betanew=breg(zmat[,1],zmat[,2:(k+1)],betabar,A)
   #
   # draw sigmai given w and beta
   #
   epsilon=matrix((wnew-X%*%betanew),nrow=pm1)
   S=crossprod(t(epsilon))
   W=rwishart(nu+n,chol2inv(chol(V+S)))
   C=W$C
   #
   #       print time to completion and draw # every 100th draw
   #
   if(rep%%100 == 0)
     {ctime=proc.time()[3]
      timetoend=((ctime-itime)/rep)*(R+1-rep)
      cat(" ",rep," (",round(timetoend/60,1),")",fill=TRUE)
      fsh()}
   #
   #       save every keepth draw
   #
   if(rep%%keep ==0)
      {mkeep=rep/keep
      betadraw[mkeep,]=betanew
      sigmadraw[mkeep,]=as.vector(W$IW)}
   wold=wnew
   betaold=betanew
   }
ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')

attributes(betadraw)$class=c("bayesm.mat","mcmc")
attributes(betadraw)$mcpar=c(1,R,keep)
attributes(sigmadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(sigmadraw)$mcpar=c(1,R,keep)
list(betadraw=betadraw,sigmadraw=sigmadraw)
}
