rhierNegbinRw = 
function(Data, Prior, Mcmc) {

#   Revision History
#	  Sridhar Narayanan - 05/2005
#         P. Rossi 6/05
#         fixed error with nobs not specified and changed llnegbinFract 9/05
#         3/07 added classes
#
#   Model
#       (y_i|lambda_i,alpha) ~ Negative Binomial(Mean = lambda_i, Overdispersion par = alpha)
#
#       ln(lambda_i) =  X_i * beta_i
#
#       beta_i = Delta'*z_i + nu_i
#               nu_i~N(0,Vbeta)
#
#   Priors
#       vec(Delta|Vbeta) ~ N(vec(Deltabar), Vbeta (x) (Adelta^-1))
#       Vbeta ~ Inv Wishart(nu, V)
#       alpha ~ Gamma(a,b) where mean = a/b and variance = a/(b^2)
#
#   Arguments
#       Data = list of regdata,Z 
#           regdata is a list of lists each list with members y, X
#              e.g. regdata[[i]]=list(y=y,X=X)
#              X has nvar columns including a first column of ones
#              Z is nreg=length(regdata) x nz with a first column of ones
#
#       Prior - list containing the prior parameters
#           Deltabar, Adelta - mean of Delta prior, inverse of variance covariance of Delta prior
#           nu, V - parameters of Vbeta prior
#           a, b - parameters of alpha prior
#
#       Mcmc - list containing
#           R is number of draws
#           keep is thinning parameter (def = 1)
#           s_beta - scaling parameter for beta RW (def = 2.93/sqrt(nvar))
#           s_alpha - scaling parameter for alpha RW (def = 2.93)
#           c - fractional weighting parameter (def = 2)
#           Vbeta0, Delta0 - initial guesses for parameters, if not supplied default values are used
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
     out=.Internal(dnbinom(y,alpha,prob,TRUE))
     return(sum(out))
}

llnegbinFract = 
function(par,X,y,Xpooled, ypooled, power, nvar,lnalpha)  {
# Computes the fractional log-likelihood at the unit level
#    = l_i * l_bar^power
    theta = c(par,lnalpha)
    llnegbin(theta,X,y,nvar) + power*llnegbin(theta,Xpooled,ypooled, nvar) 
}

lpostbetai = 
function(beta, alpha, X, y, Delta, Z, Vbetainv) {
# Computes the unnormalized log posterior for beta at the unit level
    lambda = exp(X %*% as.vector(beta))
    p = alpha/(alpha + lambda)
    residual = as.vector(beta - as.vector(Z%*%Delta))
    sum(alpha * log(p) + y * log(1-p)) - 0.5*( t(residual)%*%Vbetainv%*%residual)
}


lpostalpha = 
function(alpha, beta, regdata, ypooled, a, b, nreg) {
# Computes the unnormalized log posterior for alpha
    Xbeta=NULL
    for (i in 1:nreg) {Xbeta = rbind(Xbeta,regdata[[i]]$X%*%beta[i,]) }
    sum(log(dnbinom(ypooled,size=alpha,mu=exp(Xbeta)))) + (a-1)*log(alpha) - b* alpha
}


#
# Error Checking
#
pandterm=function(message) {stop(message,call.=FALSE)}
if(missing(Data)) {pandterm("Requires Data argument -- list of regdata and (possibly) Z")}

if(is.null(Data$regdata)) {
    pandterm("Requires Data element regdata -- list of data for each unit : y and X")
}
regdata=Data$regdata
nreg = length(regdata)

if (is.null(Data$Z)) {
    cat("Z not specified - using a column of ones instead", fill = TRUE)
    Z = matrix(rep(1,nreg),ncol=1)
}
else {
    if (nrow(Data$Z) != nreg) {
        pandterm(paste("Nrow(Z) ", nrow(Z), "ne number units ",nreg))
    }
    else {
        Z = Data$Z
    }
}
nz = ncol(Z)

dimfun = function(l) {
    c(length(l$y),dim(l$X))
}
dims=sapply(regdata,dimfun)
dims = t(dims)
nvar = quantile(dims[,3],prob=0.5)
for (i in 1:nreg) {
        if (dims[i, 1] != dims[i, 2] || dims[i, 3] != nvar) {
            pandterm(paste("Bad Data dimensions for unit ", i, 
                " dims(y,X) =", dims[i, ]))
        }
}

ypooled = NULL
Xpooled = NULL
for (i in 1:nreg) {
    ypooled = c(ypooled,regdata[[i]]$y)
    Xpooled = rbind(Xpooled,regdata[[i]]$X)
}
nobs= length(ypooled)

nvar=ncol(Xpooled)
#
# check for prior elements
#
if(missing(Prior)) {
    Deltabar=matrix(rep(0,nvar*nz),nrow=nz) ; Adelta=0.01*diag(nz) ; nu=nvar+3; V=nu*diag(nvar); a=0.5; b=0.1;
}
else {
    if(is.null(Prior$Deltabar)) {Deltabar=matrix(rep(0,nvar*nz),nrow=nz)} else {Deltabar=Prior$Deltabar}
    if(is.null(Prior$Adelta)) {Adelta=0.01*diag(nz)} else {Adelta=Prior$Adelta}
    if(is.null(Prior$nu)) {nu=nvar+3} else {nu=Prior$nu}
    if(is.null(Prior$V)) {V=nu*diag(nvar)} else {V=Prior$V}
    if(is.null(Prior$a)) {a=0.5} else {a=Prior$a}
    if(is.null(Prior$b)) {b=0.1} else {b=Prior$b}
}

if(sum(dim(Deltabar) == c(nz,nvar)) != 2) pandterm("Deltabar is of incorrect dimension")
if(sum(dim(Adelta)==c(nz,nz)) != 2) pandterm("Adelta is of incorrect dimension")
if(nu < nvar) pandterm("invalid nu value")
if(sum(dim(V)==c(nvar,nvar)) != 2) pandterm("V is of incorrect dimension")
if((length(a) != 1) | (a <=0)) pandterm("a should be a positive number")
if((length(b) != 1) | (b <=0)) pandterm("b should be a positive number")

#
# check for Mcmc 
#
if(missing(Mcmc)) pandterm("Requires Mcmc argument -- at least R")
if(is.null(Mcmc$R)) {pandterm("Requires element R of Mcmc")} else {R=Mcmc$R}
if(is.null(Mcmc$Vbeta0)) {Vbeta0=diag(nvar)} else {Vbeta0=Mcmc$Vbeta0}
if(sum(dim(Vbeta0) == c(nvar,nvar)) !=2) pandterm("Vbeta0 is not of dimension nvar")
if(is.null(Mcmc$Delta0)) {Delta0=matrix(rep(0,nz*nvar),nrow=nz)} else {Delta0=Mcmc$Delta0}
if(sum(dim(Delta0) == c(nz,nvar)) !=2) pandterm("Delta0 is not of dimension nvar by nz")
if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
if(is.null(Mcmc$s_alpha)) {cat("Using default s_alpha = 2.93",fill=TRUE); s_alpha=2.93} 
    else {s_alpha= Mcmc$s_alpha }
if(is.null(Mcmc$s_beta)) {cat("Using default s_beta = 2.93/sqrt(nvar)",fill=TRUE); s_beta=2.93/sqrt(nvar)} 
    else {s_beta=Mcmc$s_beta }
if(is.null(Mcmc$c)) {cat("Using default c = 2"); c=2} 
    else {c = Mcmc$c}

#out = rhierNegbinRw(Data, Prior, Mcmc)
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting Random Walk Metropolis Sampler for Hierarchical Negative Binomial Regression",fill=TRUE)
cat("  ",nobs," obs; ",nvar," covariates (including the intercept); ",fill=TRUE)
cat("  ",nz," individual characteristics (including the intercept) ",fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parameters:",fill=TRUE)
cat("Deltabar",fill=TRUE)
print(Deltabar)
cat("Adelta",fill=TRUE)
print(Adelta)
cat("nu",fill=TRUE)
print(nu)
cat("V",fill=TRUE)
print(V)
cat("a",fill=TRUE)
print(a)
cat("b",fill=TRUE)
print(b)
cat(" ",fill=TRUE)
cat("MCMC Parameters:",fill=TRUE)
cat(R," reps; keeping every ",keep,"th draw",fill=TRUE)
cat("s_alpha = ",s_alpha,fill=TRUE)
cat("s_beta = ",s_beta,fill=TRUE)
cat("Fractional Scaling c = ",c,fill=TRUE)
cat(" ",fill=TRUE)

par = rep(0,(nvar+1))
cat("initializing Metropolis candidate densities for ",nreg,"units ...",fill=TRUE)
fsh()
mle = optim(par,llnegbin, X=Xpooled, y=ypooled, nvar=nvar, 
      method="L-BFGS-B", upper=c(Inf,Inf,Inf,log(100000000)), hessian=TRUE, control=list(fnscale=-1))
fsh()
beta_mle=mle$par[1:nvar]
alpha_mle = exp(mle$par[nvar+1])
varcovinv = -mle$hessian
Delta = Delta0
Beta = t(matrix(rep(beta_mle,nreg),ncol=nreg))
Vbetainv = solve(Vbeta0)
Vbeta = Vbeta0
alpha = alpha_mle
alphacvar = s_alpha/varcovinv[nvar+1,nvar+1]
alphacroot = sqrt(alphacvar)
#cat("beta_mle = ",beta_mle,fill=TRUE)
#cat("alpha_mle = ",alpha_mle, fill = TRUE)
#fsh()

hess_i=NULL
# Find the individual candidate hessian
for (i in 1:nreg) {
    power = length(regdata[[i]]$y)/(c*nobs)
    mle2 = optim(mle$par[1:nvar],llnegbinFract, X=regdata[[i]]$X, y=regdata[[i]]$y, Xpooled=Xpooled, 
           ypooled=ypooled, power=power, nvar=nvar, lnalpha=mle$par[nvar+1], 
           method="BFGS", hessian=TRUE, control=list(fnscale=-1, trace=0))
    if (mle2$convergence==0)
        hess_i[[i]] = list(hess=-mle2$hessian)
    else
        hess_i[[i]] = diag(rep(0,nvar))
   if(i%%50 ==0) cat("  completed unit #",i,fill=TRUE)	
   fsh()
}

oldlpostbeta = rep(0,nreg)
nacceptbeta = 0
nacceptalpha = 0
clpostbeta = rep(0,nreg)

Betadraw = array(double((floor(R/keep)) * nreg * nvar), dim = c(nreg, 
        nvar, floor(R/keep)))

alphadraw = rep(0,floor(R/keep))
llike = rep(0,floor(R/keep))
Vbetadraw=matrix(double(floor(R/keep)*(nvar*nvar)),ncol=(nvar*nvar))
Deltadraw=matrix(double(floor(R/keep)*(nvar*nz)),ncol=(nvar*nz))

#
#	start main iteration loop
#
itime=proc.time()[3]
cat(" ",fill=TRUE)
cat("MCMC Iteration (est time to end - min) ",fill=TRUE)
fsh()

for (r in 1:R) 
{
#   Draw betai
    for (i in 1:nreg) {
        betacvar = s_beta*solve(hess_i[[i]]$hess + Vbetainv)
        betaroot = t(chol(betacvar))
        betac = as.vector(Beta[i,]) + betaroot%*%rnorm(nvar)

        oldlpostbeta[i] = lpostbetai(as.vector(Beta[i,]), alpha, regdata[[i]]$X, regdata[[i]]$y, Delta, Z[i,],Vbetainv)
        clpostbeta[i] = lpostbetai(betac, alpha, regdata[[i]]$X, regdata[[i]]$y, Delta, Z[i,],Vbetainv)
        
        ldiff=clpostbeta[i]-oldlpostbeta[i]
        acc=min(1,exp(ldiff))
        if(acc < 1) {unif=runif(1)} else {unif=0}

        if (unif <= acc) {
            Beta[i,]=betac
            nacceptbeta=nacceptbeta+1
        }
    }

#   Draw alpha
    logalphac = rnorm(1,mean=log(alpha), sd=alphacroot)
    oldlpostalpha = lpostalpha(alpha, Beta, regdata, ypooled,  a, b, nreg)
    clpostalpha = lpostalpha(exp(logalphac), Beta, regdata, ypooled, a, b, nreg)
    ldiff=clpostalpha-oldlpostalpha
    acc=min(1,exp(ldiff))
    if(acc < 1) {unif=runif(1)} else {unif=0}
    if (unif <= acc) {
        alpha=exp(logalphac)
        nacceptalpha=nacceptalpha+1
    }

#   Draw Vbeta and Delta using rmultireg (bayesm function)
    temp = rmultireg(Beta,Z,Deltabar,Adelta,nu,V)
    Vbeta = matrix(temp$Sigma,nrow=nvar)
    Vbetainv = solve(Vbeta)
    Delta = temp$B


  if(r%%100 == 0)
    {ctime=proc.time()[3]
    timetoend=((ctime-itime)/r)*(R-r)
    cat(" ",r," (",round(timetoend/60,1),")",fill=TRUE)
    fsh()}

  if(r%%keep == 0) {
    mkeep=r/keep
    Betadraw[, ,mkeep]=Beta
    alphadraw[mkeep] = alpha
    Vbetadraw[mkeep,] = as.vector(Vbeta)
    Deltadraw[mkeep,] = as.vector(Delta)
    ll=0.0
    for (i in 1:nreg) {ll=ll+llnegbin(c(Beta[i,],alpha),regdata[[i]]$X,regdata[[i]]$y,nvar)}
    llike[r]=ll
  }
}
ctime = proc.time()[3]

attributes(alphadraw)$class=c("bayesm.mat","mcmc")
attributes(alphadraw)$mcpar=c(1,R,keep)
attributes(Deltadraw)$class=c("bayesm.mat","mcmc")
attributes(Deltadraw)$mcpar=c(1,R,keep)
attributes(Vbetadraw)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(Vbetadraw)$mcpar=c(1,R,keep)
attributes(Betadraw)$class=c("bayesm.hcoef")
    
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')
return(list(llike=llike,Betadraw=Betadraw,alphadraw=alphadraw, Vbetadraw=Vbetadraw, Deltadraw=Deltadraw,
     acceptrbeta=nacceptbeta/(R*nreg)*100,acceptralpha=nacceptalpha/R*100))
}
