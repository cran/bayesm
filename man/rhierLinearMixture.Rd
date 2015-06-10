\name{rhierLinearMixture}
\alias{rhierLinearMixture}
\concept{bayes}
\concept{MCMC}
\concept{Gibbs Sampling}
\concept{mixture of normals}
\concept{normal mixture}
\concept{heterogeneity}
\concept{regresssion}
\concept{hierarchical models}
\concept{linear model}

\title{ Gibbs Sampler for Hierarchical Linear Model }
\description{
  \code{rhierLinearMixture} implements a Gibbs Sampler for hierarchical linear models with a mixture of normals prior.
}
\usage{
rhierLinearMixture(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(regdata,Z) (Z optional). }
  \item{Prior}{ list(deltabar,Ad,mubar,Amu,nu,V,nu.e,ssq,ncomp)  (all but ncomp are optional).}
  \item{Mcmc}{ list(R,keep,nprint) (R required).}
}
\details{
  Model: length(regdata) regression equations. \cr
        \eqn{y_i = X_i\beta_i + e_i}. \eqn{e_i} \eqn{\sim}{~} \eqn{N(0,\tau_i)}.  \code{nvar} is the number of X vars in each equation. 

 Priors:\cr
        \eqn{\tau_i} \eqn{\sim}{~} \eqn{nu.e*ssq_i/\chi^2_{nu.e}}.  \eqn{\tau_i} is the variance of \eqn{e_i}.\cr
        \eqn{B = Z\Delta + U} or \cr
        \eqn{\beta_i = \Delta' Z[i,]' + u_i}. \cr
        \eqn{\Delta} is an nz x nvar array. \cr


        \eqn{u_i} \eqn{\sim}{~} \eqn{N(\mu_{ind},\Sigma_{ind})}\cr
        \eqn{ind} \eqn{\sim}{~} \eqn{multinomial(pvec)} \cr

        \eqn{pvec} \eqn{\sim}{~} \eqn{dirichlet(a)}\cr
        \eqn{delta= vec(\Delta)} \eqn{\sim}{~} \eqn{N(deltabar,A_d^{-1})}\cr
        \eqn{\mu_j} \eqn{\sim}{~} \eqn{N(mubar,\Sigma_j (x) Amu^{-1})}\cr
        \eqn{\Sigma_j} \eqn{\sim}{~} \eqn{IW(nu,V)} \cr


  List arguments contain:
  \itemize{
    \item{\code{regdata}}{ list of lists with X,y matrices for each of length(regdata) regressions}
    \item{\code{regdata[[i]]$X}}{ X matrix for equation i }
    \item{\code{regdata[[i]]$y}}{ y vector for equation i }
    \item{\code{deltabar}}{nz*nvar vector of prior means (def: 0)}
    \item{\code{Ad}}{ prior prec matrix for vec(Delta) (def: .01I)}
    \item{\code{mubar}}{ nvar x 1 prior mean vector for normal comp mean (def: 0)}
    \item{\code{Amu}}{ prior precision for normal comp mean (def: .01I)}
    \item{\code{nu}}{ d.f. parm for IW prior on norm comp Sigma (def: nvar+3)}
    \item{\code{V}}{ pds location parm for IW prior on norm comp Sigma (def: nuI)}
    \item{\code{nu.e}}{ d.f. parm for regression error variance prior (def: 3)}
    \item{\code{ssq}}{ scale parm for regression error var prior (def: var(\eqn{y_i}))}
    \item{\code{a}}{ Dirichlet prior parameter (def: 5)}
    \item{\code{ncomp}}{ number of components used in normal mixture }
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
   }
}
\value{
  a list containing
  \item{taudraw}{R/keep x nreg array of error variance draws}
  \item{betadraw}{nreg x nvar x R/keep array of individual regression coef draws}
  \item{Deltadraw}{R/keep x nz x nvar array of Deltadraws}
  \item{nmix}{list of three elements, (probdraw, NULL, compdraw)}
}
\note{
  More on \code{probdraw} component of nmix return value list: \cr
  this is an R/keep by ncomp array of draws of mixture component probs (pvec)\cr
  More on \code{compdraw} component of nmix return value list: 
  \describe{
  \item{compdraw[[i]]}{the ith draw of components for mixtures}
  \item{compdraw[[i]][[j]]}{ith draw of the jth normal mixture comp}
  \item{compdraw[[i]][[j]][[1]]}{ith draw of jth normal mixture comp mean vector}
  \item{compdraw[[i]][[j]][[2]]}{ith draw of jth normal mixture cov parm (rooti)}
  }

  Note: Z should \strong{not} include an intercept and should be centered for ease of interpretation. The mean of each of the \code{nreg} \eqn{\beta} s is the mean of the normal mixture.  Use \code{summary()} to compute this mean from the \code{compdraw} output.  
  
  
  Be careful in assessing the prior parameter, Amu.  .01 can be too small for some applications. See 
  Rossi et al, chapter 5 for full discussion.\cr

} 
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 3. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\seealso{ \code{\link{rhierLinearModel}} }
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(66)
nreg=300; nobs=500; nvar=3; nz=2

Z=matrix(runif(nreg*nz),ncol=nz) 
Z=t(t(Z)-apply(Z,2,mean))
Delta=matrix(c(1,-1,2,0,1,0),ncol=nz)
tau0=.1
iota=c(rep(1,nobs))

## create arguments for rmixture

tcomps=NULL
a=matrix(c(1,0,0,0.5773503,1.1547005,0,-0.4082483,0.4082483,1.2247449),ncol=3)
tcomps[[1]]=list(mu=c(0,-1,-2),rooti=a) 
tcomps[[2]]=list(mu=c(0,-1,-2)*2,rooti=a)
tcomps[[3]]=list(mu=c(0,-1,-2)*4,rooti=a)
tpvec=c(.4,.2,.4)                               

regdata=NULL						  # simulated data with Z
betas=matrix(double(nreg*nvar),ncol=nvar)
tind=double(nreg)

for (reg in 1:nreg) {
tempout=rmixture(1,tpvec,tcomps)
betas[reg,]=Delta\%*\%Z[reg,]+as.vector(tempout$x)
tind[reg]=tempout$z
X=cbind(iota,matrix(runif(nobs*(nvar-1)),ncol=(nvar-1)))
tau=tau0*runif(1,min=0.5,max=1)
y=X\%*\%betas[reg,]+sqrt(tau)*rnorm(nobs)
regdata[[reg]]=list(y=y,X=X,beta=betas[reg,],tau=tau)
}

## run rhierLinearMixture

Data1=list(regdata=regdata,Z=Z)
Prior1=list(ncomp=3)
Mcmc1=list(R=R,keep=1)

out1=rhierLinearMixture(Data=Data1,Prior=Prior1,Mcmc=Mcmc1)

cat("Summary of Delta draws",fill=TRUE)
summary(out1$Deltadraw,tvalues=as.vector(Delta))
cat("Summary of Normal Mixture Distribution",fill=TRUE)
summary(out1$nmix)

if(0){
## plotting examples 
plot(out1$betadraw)
plot(out1$nmix)
plot(out1$Deltadraw)
}

}
\keyword{ regression }
