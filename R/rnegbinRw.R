rnegbinRw = 
function(Data, Prior, Mcmc) {

#   Revision History
#	  Sridhar Narayanan - 05/2005
#         P. Rossi 6/05
#         3/07 added classes
#
#   Model
#       (y|lambda,alpha) ~ Negative Binomial(Mean = lambda, Overdispersion par = alpha)
#
#       ln(lambda) =  X * beta
#               
#   Priors
#       beta ~ N(betabar, A^-1)
#       alpha ~ Gamma(a,b) where mean = a/b and variance = a/(b^2)
#
#   Arguments
#       Data = list of y, X
#              e.g. regdata[[i]]=list(y=y,X=X)
#              X has nvar columns including a first column of ones
#
#       Prior - list containing the prior parameters
#           betabar, A - mean of beta prior, inverse of variance covariance of beta prior
#           a, b - parameters of alpha prior
#
#       Mcmc - list containing
#           R is number of draws
#           keep is thinning parameter (def = 1)
#           s_beta - scaling parameter for beta RW (def = 2.93/sqrt(nvar))
#           s_alpha - scaling parameter for alpha RW (def = 2.93)
#           beta0 - initial guesses for parameters, if not supplied default values are used
#


#
# Definitions of functions used within rhierNegbinRw
#
llnegbin = 
function(par,X,y, nvar) {
# Computes the log-likelihood
    beta = par[1:nvar]
    alpha = exp(par[nvar+1])+1.0e-50
    mean=exp(X%*%beta)
    prob=alpha/(alpha+mean)
    prob=ifelse(prob<1.0e-100,1.0e-100,prob)
     out=dnbinom(y,size=alpha,prob=prob,log=TRUE)
     return(sum(out))
}

lpostbetai = 
function(beta, alpha, X, y, betabar, A) {
# Computes the unnormalized log posterior for beta
    lambda = exp(X %*% beta)
    p = alpha/(alpha + lambda)
    residual = as.vector(beta - betabar)
    sum(alpha * log(p) + y * log(1-p)) - 0.5*( t(residual)%*%A%*%residual)
}


lpostalpha = 
function(alpha, beta, X,y, a, b) {
# Computes the unnormalized log posterior for alpha
    sum(log(dnbinom(y,size=alpha,mu=exp(X%*%beta)))) + (a-1)*log(alpha) - b* alpha
}


#
# Error Checking
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of X and y")}
if(is.null(Data$X)) {pandterm("Requires Data element X")} else {X=Data$X}
if(is.null(Data$y)) {pandterm("Requires Data element y")} else {y=Data$y}
nvar = ncol(X)

if (length(y) != nrow(X)) {pandterm("Mismatch in the number of observations in X and y")}
nobs=length(y)

#
# check for prior elements
#
if(missing(Prior)) {
    betabar=rep(0,nvar); A=0.01*diag(nvar) ;  a=0.5; b=0.1;
}
else {
    if(is.null(Prior$betabar)) {betabar=rep(0,nvar)} else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=0.01*diag(nvar)} else {A=Prior$A}
    if(is.null(Prior$a)) {a=0.5} else {a=Prior$a}
    if(is.null(Prior$b)) {b=0.1} else {b=Prior$b}
}

if(length(betabar) != nvar) pandterm("betabar is of incorrect dimension")
if(sum(dim(A)==c(nvar,nvar)) != 2) pandterm("A is of incorrect dimension")
if((length(a) != 1) | (a <=0)) pandterm("a should be a positive number")
if((length(b) != 1) | (b <=0)) pandterm("b should be a positive number")

#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$beta0)) {beta0=rep(0,nvar)} else {beta0=Mcmc$beta0}
if(length(beta0) !=nvar) pandterm("beta0 is not of dimension nvar")
if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
if(is.null(Mcmc$s_alpha)) {cat("Using default s_alpha = 2.93",fill=TRUE); s_alpha=2.93} 
    else {s_alpha = Mcmc$s_alpha} 
if(is.null(Mcmc$s_beta)) {cat("Using default s_beta = 2.93/sqrt(nvar)",fill=TRUE); s_beta=2.93/sqrt(nvar)} 
    else {s_beta = Mcmc$s_beta}

#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting Random Walk Metropolis Sampler for Negative Binomial Regression",fill=TRUE)
cat("  ",nobs," obs; ",nvar," covariates (including intercept); ",fill=TRUE)
cat("Prior Parameters:",fill=TRUE)
cat("betabar",fill=TRUE)
print(betabar)
cat("A",fill=TRUE)
print(A)
cat("a",fill=TRUE)
print(a)
cat("b",fill=TRUE)
print(b)
cat(" ",fill=TRUE)
cat("MCMC Parms: ",fill=TRUE)
cat("  ",R," reps; keeping every ",keep,"th draw",fill=TRUE)
cat("s_alpha = ",s_alpha,fill=TRUE)
cat("s_beta = ",s_beta,fill=TRUE)
cat(" ",fill=TRUE)

par = rep(0,(nvar+1))
cat(" Initializing RW Increment Covariance Matrix...",fill=TRUE)
fsh()
mle = optim(par,llnegbin, X=X, y=y, nvar=nvar, method="L-BFGS-B", upper=c(Inf,Inf,Inf,log(100000000)), hessian=TRUE, control=list(fnscale=-1))
fsh()
beta_mle=mle$par[1:nvar]
alpha_mle = exp(mle$par[nvar+1])
varcovinv = -mle$hessian
beta = beta0
betacvar = s_beta*solve(varcovinv[1:nvar,1:nvar])
betaroot = t(chol(betacvar))
alpha = alpha_mle
alphacvar = s_alpha/varcovinv[nvar+1,nvar+1]
alphacroot = sqrt(alphacvar)
cat("beta_mle = ",beta_mle,fill=TRUE)
cat("alpha_mle = ",alpha_mle, fill = TRUE)
fsh()

oldlpostbeta = 0
nacceptbeta = 0
nacceptalpha = 0
clpostbeta = 0

alphadraw = rep(0,floor(R/keep))
betadraw=matrix(double(floor(R/keep)*(nvar)),ncol=nvar)
llike=rep(0,floor(R/keep))

#
#	start main iteration loop
#
itime=proc.time()[3]
cat(" ",fill=TRUE)
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
fsh()

for (r in 1:R) 
{
#   Draw beta
    betac = beta + betaroot%*%rnorm(nvar)
    oldlpostbeta = lpostbetai(beta, alpha, X, y, betabar,A)
    clpostbeta = lpostbetai(betac, alpha, X, y, betabar,A)
        
    ldiff=clpostbeta-oldlpostbeta
    acc=min(1,exp(ldiff))
    if(acc < 1) {unif=runif(1)} else {unif=0}
    if (unif <= acc) {
        beta=betac
        nacceptbeta=nacceptbeta+1
    }
 

#   Draw alpha
    logalphac = rnorm(1,mean=log(alpha), sd=alphacroot)
    oldlpostalpha = lpostalpha(alpha, beta, X, y,  a, b)
    clpostalpha = lpostalpha(exp(logalphac), beta, X, y, a, b)
    ldiff=clpostalpha-oldlpostalpha
    acc=min(1,exp(ldiff))
    if(acc < 1) {unif=runif(1)} else {unif=0}
    if (unif <= acc) {
        alpha=exp(logalphac)
        nacceptalpha=nacceptalpha+1
    }

if(r%%100 == 0)
    {ctime=proc.time()[3]
    timetoend=((ctime-itime)/r)*(R-r)
    cat(" ",r," (",round(timetoend/60,1),")",fill=TRUE)
    fsh()}

if(r%%keep == 0) {
    mkeep=r/keep
    betadraw[mkeep,]=beta
    alphadraw[mkeep] = alpha
    llike[mkeep]=llnegbin(c(beta,alpha),X,y,nvar)
    }
}

ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')


attributes(betadraw)$class=c("bayesm.mat","mcmc")
attributes(betadraw)$mcpar=c(1,R,keep)
attributes(alphadraw)$class=c("bayesm.mat","mcmc")
attributes(alphadraw)$mcpar=c(1,R,keep)
return(list(llike=llike,betadraw=betadraw,alphadraw=alphadraw,
     acceptrbeta=nacceptbeta/R*100,acceptralpha=nacceptalpha/R*100))
}
