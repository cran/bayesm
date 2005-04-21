rmvpGibbs=
function(Data,Prior,Mcmc) 
{
#
# Revision History:
#   modified by rossi 12/18/04 to include error checking
#
# purpose:  Gibbs MVP model with full covariance matrix
#
# Arguments:
#   Data contains 
#      p the number of alternatives (could be time or could be from pick j of p survey)
#      y -- a vector of length n*p of indicators (1 if "chosen" if not)
#      X -- np x k matrix of covariates (including intercepts)
#                 each X_i is p x nvar
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
#    w_i = X_ibeta + e    e~N(0,Sigma)     note w_i,e are p x 1
#    y_ij = 1 if w_ij > 0 else y_ij = 0  
#  
#  priors:
#    beta ~ N(betabar,A^-1) in prior
#    Sigma ~ IW(nu,V)
#
#  Check arguments
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of p, y, X")}
  if(is.null(Data$p)) {pandterm("Requires Data element p -- number of binary indicators")}
  p=Data$p
  if(is.null(Data$y)) {pandterm("Requires Data element y -- values of binary indicators")}
  y=Data$y
  if(is.null(Data$X)) {pandterm("Requires Data element X -- matrix of covariates")}
  X=Data$X
#
# check data for validity
#
levely=as.numeric(levels(as.factor(y)))
  bady=FALSE
  for (i in 0:1) 
  { if(levely[i+1] != i) {bady=TRUE} }
cat("Table of y values",fill=TRUE)
print(table(y))
if (bady) {pandterm("Invalid y")}
if (length(y)%%p !=0) {pandterm("length of y is not a multiple of p")}
n=length(y)/p
k=ncol(X)
if(nrow(X) != (n*p)) {pandterm(paste("X has ",nrow(X)," rows; must be = p*n"))}
#
# check for prior elements
#
if(missing(Prior)) 
  { betabar=rep(0,k) ; A=.01*diag(k) ; nu=p+3; V=nu*diag(p)}
else 
  {if(is.null(Prior$betabar)) {betabar=rep(0,k)} else {betabar=Prior$betabar}
   if(is.null(Prior$A)) {A=.01*diag(k)} else {A=Prior$A}
   if(is.null(Prior$nu)) {nu=p+3} else {nu=Prior$nu}
   if(is.null(Prior$V)) {V=nu*diag(p)} else {V=Prior$V}}
if(length(betabar) != k) pandterm("length betabar ne k")
if(sum(dim(A)==c(k,k)) != 2) pandterm("A is of incorrect dimension")
if(nu < 1) pandterm("invalid nu value")
if(sum(dim(V)==c(p,p)) != 2) pandterm("V is of incorrect dimension")
#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R must be included")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$beta0)) {beta0=rep(0,k)} else {beta0=Mcmc$beta0}
if(is.null(Mcmc$sigma0)) {sigma0=diag(p)} else {sigma0=Mcmc$sigma0}
if(length(beta0) != k) pandterm("beta0 is not of length k")
if(sum(dim(sigma0) == c(p,p)) != 2) pandterm("sigma0 is of incorrect dimension")
if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting Gibbs Sampler for MVP",fill=TRUE)
cat("  ",n," obs of ",p," binary indicators; ",k," indep vars (including intercepts)",fill=TRUE)
cat("  ",R," reps; keeping every ",keep,"th draw",fill=TRUE)
cat(" ",fill=TRUE)
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
sigmadraw=matrix(double(floor(R/keep)*p*p),ncol=p*p)
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
drawwMvpC=function(w,mu,y,sigi) {
	p=ncol(sigi)
      .C("draww_mvp",w=as.double(w),as.double(mu),as.double(sigi),
        as.integer(length(w)/p),as.integer(p),as.integer(y))$w}

drawwMvp=
function(w,X,y,beta,sigmai){
#
#   draw latent vector
#
#  	w is n x (p-1) vector
#       X ix n(p-1) x k  matrix
#       y is n x (p-1) vector of binary (0,1) outcomes 
#       beta is k x 1 vector
#       sigmai is (p-1) x (p-1) 
#

Xbeta=as.vector(X%*%beta)
drawwMvpC(w,Xbeta,y,sigmai)
}

itime=proc.time()[3]
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
for (rep in 1:R) 
   {
   #
   # draw w given beta(rep-1),sigma(rep-1)
   #
   sigmai=crossprod(C)
   wnew=drawwMvp(wold,X,y,betaold,sigmai)
   #
   # draw beta given w(rep) and sigma(rep-1)
   #
   #  note:  if Sigma^-1 (G) = C'C then Var(Ce)=CSigmaC' = I
   #  first, transform w_i = X_ibeta + e_i by premultiply by C
   #
   zmat=matrix(cbind(wnew,X),nrow=p)
   zmat=C%*%zmat
   zmat=matrix(zmat,nrow=nrow(X))
   betanew=breg(zmat[,1],zmat[,2:(k+1)],betabar,A)
   #
   # draw sigmai given w and beta
   #
   epsilon=matrix((wnew-X%*%betanew),nrow=p)
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

return(list(betadraw=betadraw,sigmadraw=sigmadraw))
}
