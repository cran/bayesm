\name{rhierMnlRwMixture}
\alias{rhierMnlRwMixture}
\concept{bayes}
\concept{MCMC}
\concept{Multinomial Logit}
\concept{mixture of normals}
\concept{normal mixture}
\concept{heterogeneity}
\concept{hierarchical models}

\title{ MCMC Algorithm for Hierarchical Multinomial Logit with Mixture of Normals Heterogeneity}
\description{
  \code{rhierMnlRwMixture} is a MCMC algorithm for a hierarchical multinomial logit with a mixture of normals 
  heterogeneity distribution.  This is a hybrid Gibbs Sampler with a RW Metropolis step for the MNL 
  coefficients for each panel unit.
}
\usage{
rhierMnlRwMixture(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(p,lgtdata,Z) ( Z is optional) }
  \item{Prior}{ list(a,deltabar,Ad,mubar,Amu,nu,V,a,ncomp) (all but ncomp are optional)}
  \item{Mcmc}{ list(s,w,R,keep,nprint) (R required)}
}
\details{
  Model: \cr
  \eqn{y_i} \eqn{\sim}{~} \eqn{MNL(X_i,\beta_i)}.  i=1,\ldots, length(lgtdata). \eqn{\beta_i} is nvar x 1.

  \eqn{\beta_i}= Z\eqn{\Delta}[i,] + \eqn{u_i}. \cr
  Note:  Z\eqn{\Delta} is the matrix Z * \eqn{\Delta}; [i,] refers to ith row of this product.\cr
  Delta is an nz x nvar array. 

  \eqn{u_i} \eqn{\sim}{~} \eqn{N(\mu_{ind},\Sigma_{ind})}. \eqn{ind} \eqn{\sim}{~} multinomial(pvec). \cr

  Priors: \cr
  \eqn{pvec} \eqn{\sim}{~} dirichlet (a)\cr
  \eqn{delta= vec(\Delta)} \eqn{\sim}{~} \eqn{N(deltabar,A_d^{-1})}\cr
  \eqn{\mu_j} \eqn{\sim}{~} \eqn{N(mubar,\Sigma_j (x) Amu^{-1})}\cr
  \eqn{\Sigma_j} \eqn{\sim}{~} IW(nu,V) \cr

  Lists contain:
  \itemize{
    \item{\code{p}}{ p is number of choice alternatives}
    \item{\code{lgtdata}}{list of lists with each cross-section unit MNL data}
    \item{\code{lgtdata[[i]]$y}}{ \eqn{n_i} vector of multinomial outcomes (1,\ldots,m)}
    \item{\code{lgtdata[[i]]$X}}{ \eqn{n_i}*p by nvar design matrix for ith unit}
    \item{\code{a}}{vector of length ncomp of Dirichlet prior parms (def: rep(5,ncomp))}
    \item{\code{deltabar}}{nz*nvar vector of prior means (def: 0)}
    \item{\code{Ad}}{ prior prec matrix for vec(D) (def: .01I)}
    \item{\code{mubar}}{ nvar x 1 prior mean vector for normal comp mean (def: 0)}
    \item{\code{Amu}}{ prior precision for normal comp mean (def: .01I)}
    \item{\code{nu}}{ d.f. parm for IW prior on norm comp Sigma (def: nvar+3)}
    \item{\code{V}}{ pds location parm for IW prior on norm comp Sigma (def: nuI)}
    \item{\code{a}}{ Dirichlet prior parameter (def: 5)}
    \item{\code{ncomp}}{ number of components used in normal mixture }
    \item{\code{s}}{ scaling parm for RW Metropolis (def: 2.93/sqrt(nvar))}
    \item{\code{w}}{ fractional likelihood weighting parm (def: .1)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
  }
}
\value{
  a list containing:
  \item{Deltadraw}{R/keep  x nz*nvar matrix of draws of Delta, first row is initial value}
  \item{betadraw}{ nlgt x nvar x R/keep array of draws of betas}
  \item{nmix}{ list of 3 components, probdraw, NULL, compdraw }
  \item{loglike}{ log-likelihood for each kept draw (length R/keep)}
}
\note{
  More on \code{probdraw} component of nmix list:\cr 
  R/keep x ncomp matrix of draws of probs of mixture components (pvec)  \cr
  More on \code{compdraw} component of return value list: \cr
  \itemize{
  \item{compdraw[[i]]}{ the ith draw of components for mixtures}
  \item{compdraw[[i]][[j]]}{ ith draw of the jth normal mixture comp}
  \item{compdraw[[i]][[j]][[1]]}{ ith draw of jth normal mixture comp mean vector}
  \item{compdraw[[i]][[j]][[2]]}{ ith draw of jth normal mixture cov parm (rooti) }
  }

  Note: Z should \strong{not} include an intercept and is centered for ease of interpretation. The mean of each of the \code{nlgt} \eqn{\beta} s is the mean of the normal mixture.  Use \code{summary()} to compute this mean from the \code{compdraw} output.\cr
  
  Be careful in assessing prior parameter, Amu.  .01 is too small for many applications. See 
  Rossi et al, chapter 5 for full discussion.\cr

  Note: as of version 2.0-2 of \code{bayesm}, the fractional weight parameter has been changed
  to a weight between 0 and 1.  w is the fractional weight on the normalized pooled likelihood.
  This differs from what is in Rossi et al chapter 5, i.e.

  \eqn{like_i^{(1-w)} x like_pooled^{((n_i/N)*w)}}

 
  Large R values may be required (>20,000).

} 
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 5. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
  
\seealso{ \code{\link{rmnlIndepMetrop}} }
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=10000} else {R=10}

set.seed(66)
p=3                                # num of choice alterns
ncoef=3  
nlgt=300                           # num of cross sectional units
nz=2
Z=matrix(runif(nz*nlgt),ncol=nz)
Z=t(t(Z)-apply(Z,2,mean))          # demean Z
ncomp=3                                # no of mixture components
Delta=matrix(c(1,0,1,0,1,2),ncol=2)
comps=NULL
comps[[1]]=list(mu=c(0,-1,-2),rooti=diag(rep(1,3)))
comps[[2]]=list(mu=c(0,-1,-2)*2,rooti=diag(rep(1,3)))
comps[[3]]=list(mu=c(0,-1,-2)*4,rooti=diag(rep(1,3)))
pvec=c(.4,.2,.4)

simmnlwX= function(n,X,beta) {
  ##  simulate from MNL model conditional on X matrix
  k=length(beta)
  Xbeta=X\%*\%beta
  j=nrow(Xbeta)/n
  Xbeta=matrix(Xbeta,byrow=TRUE,ncol=j)
  Prob=exp(Xbeta)
  iota=c(rep(1,j))
  denom=Prob\%*\%iota
  Prob=Prob/as.vector(denom)
  y=vector("double",n)
  ind=1:j
  for (i in 1:n) 
      {yvec=rmultinom(1,1,Prob[i,]); y[i]=ind\%*\%yvec}
  return(list(y=y,X=X,beta=beta,prob=Prob))
}

## simulate data
simlgtdata=NULL
ni=rep(50,300)
for (i in 1:nlgt) 
{  betai=Delta\%*\%Z[i,]+as.vector(rmixture(1,pvec,comps)$x)
   Xa=matrix(runif(ni[i]*p,min=-1.5,max=0),ncol=p)
   X=createX(p,na=1,nd=NULL,Xa=Xa,Xd=NULL,base=1)
   outa=simmnlwX(ni[i],X,betai)
   simlgtdata[[i]]=list(y=outa$y,X=X,beta=betai)
}

## plot betas
if(0){
## set if(1) above to produce plots
bmat=matrix(0,nlgt,ncoef)
for(i in 1:nlgt) {bmat[i,]=simlgtdata[[i]]$beta}
par(mfrow=c(ncoef,1))
for(i in 1:ncoef) hist(bmat[,i],breaks=30,col="magenta")
}

##   set parms for priors and Z
Prior1=list(ncomp=5)

keep=5
Mcmc1=list(R=R,keep=keep)
Data1=list(p=p,lgtdata=simlgtdata,Z=Z)

out=rhierMnlRwMixture(Data=Data1,Prior=Prior1,Mcmc=Mcmc1)

cat("Summary of Delta draws",fill=TRUE)
summary(out$Deltadraw,tvalues=as.vector(Delta))
cat("Summary of Normal Mixture Distribution",fill=TRUE)
summary(out$nmix)

if(0) {
## plotting examples
plot(out$betadraw)
plot(out$nmix)
}

}

\keyword{models}
