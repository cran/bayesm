\name{rhierBinLogit}
\alias{rhierBinLogit}
\concept{bayes}
\concept{MCMC}
\concept{hierarchical models}
\concept{binary logit}

\title{ MCMC Algorithm for Hierarchical Binary Logit }
\description{
  \code{rhierBinLogit} implements an MCMC algorithm for hierarchical binary logits with
  a normal heterogeneity distribution.  This is a hybrid sampler with a RW Metropolis step
  for unit-level logit parameters.

  \code{rhierBinLogit} is designed for use on choice-based conjoint data with partial profiles.
  The Design matrix is based on differences of characteristics between two alternatives. See
  Appendix A of \emph{Bayesian Statistics and Marketing} for details.
}
\usage{
rhierBinLogit(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(lgtdata,Z) (note: Z is optional) }
  \item{Prior}{ list(Deltabar,ADelta,nu,V) (note: all are optional)}
  \item{Mcmc}{ list(sbeta,R,keep) (note: all but R are optional)}
}
\details{
  Model: \cr
  \eqn{y_{hi} = 1} with \eqn{pr=exp(x_{hi}'beta_h)/(1+exp(x_{hi}'beta_h)}.  \eqn{beta_h} is nvar x 1.\cr
  h=1,\ldots,length(lgtdata) units or "respondents" for survey data.

  \eqn{beta_h}= ZDelta[h,] + \eqn{u_h}. \cr
  Note: here ZDelta refers to Z\%*\%Delta, ZDelta[h,] is hth row of this product.\cr
  Delta is an nz x nvar array. 

  \eqn{u_h} \eqn{\sim}{~} \eqn{N(0,V_{beta})}.  \cr

  Priors: \cr
  \eqn{delta= vec(Delta)} \eqn{\sim}{~} \eqn{N(vec(Deltabar),V_{beta} (x) ADelta^{-1})}\cr
  \eqn{V_{beta}} \eqn{\sim}{~} \eqn{IW(nu,V)}

  Lists contain:
  \itemize{
    \item{\code{lgtdata}}{list of lists with each cross-section unit MNL data}
    \item{\code{lgtdata[[h]]$y}}{ \eqn{n_h} vector of binary outcomes (0,1)}
    \item{\code{lgtdata[[h]]$X}}{ \eqn{n_h} by nvar design matrix for hth unit}
    \item{\code{Deltabar}}{nz x nvar matrix of prior means (def: 0)}
    \item{\code{ADelta}}{ prior prec matrix  (def: .01I)}
    \item{\code{nu}}{ d.f. parm for IW prior on norm comp Sigma (def: nvar+3)}
    \item{\code{V}}{ pds location parm for IW prior on norm comp Sigma (def: nuI)}
    \item{\code{sbeta}}{ scaling parm for RW Metropolis (def: .2)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
  }
}
\value{
  a list containing:
  \item{Deltadraw}{R/keep  x nz*nvar matrix of draws of Delta}
  \item{betadraw}{ nlgt x nvar x R/keep array of draws of betas}
  \item{Vbetadraw}{ R/keep x nvar*nvar matrix of draws of Vbeta}
  \item{llike}{R/keep vector of log-like values}
  \item{reject}{R/keep vector of reject rates over nlgt units}
}

\note{ Some experimentation with the Metropolis scaling paramter (sbeta) may be required. }

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 5. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}
\examples{
##  
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=10000} else {R=10}

set.seed(66)
nvar=5                           ## number of coefficients
nlgt=1000                        ## number of cross-sectional units
nobs=10                          ## number of observations per unit
nz=2                             ## number of regressors in mixing distribution

## set hyper-parameters
##     B=ZDelta + U  

Z=matrix(c(rep(1,nlgt),runif(nlgt,min=-1,max=1)),nrow=nlgt,ncol=nz)
Delta=matrix(c(-2,-1,0,1,2,-1,1,-.5,.5,0),nrow=nz,ncol=nvar)
iota=matrix(1,nrow=nvar,ncol=1)
Vbeta=diag(nvar)+.5*iota\%*\%t(iota)

## simulate data
lgtdata=NULL

for (i in 1:nlgt) 
{ beta=t(Delta)\%*\%Z[i,]+as.vector(t(chol(Vbeta))\%*\%rnorm(nvar))
  X=matrix(runif(nobs*nvar),nrow=nobs,ncol=nvar)
  prob=exp(X\%*\%beta)/(1+exp(X\%*\%beta)) 
  unif=runif(nobs,0,1)
  y=ifelse(unif<prob,1,0)
  lgtdata[[i]]=list(y=y,X=X,beta=beta)
}

Data=list(Dat=lgtdata,Demo=Z)
out=rhierBinLogit(Data=list(lgtdata=lgtdata,Z=Z),Mcmc=list(R=R))

cat(" Deltadraws ",fill=TRUE)
mat=apply(out$Deltadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(as.vector(Delta),mat); rownames(mat)[1]="delta"; print(mat)
cat(" Vbetadraws ",fill=TRUE)
mat=apply(out$Vbetadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(as.vector(Vbeta),mat); rownames(mat)[1]="Vbeta"; print(mat)

if(0){
td=as.vector(Delta)
par(mfrow=c(2,2))
matplot(out$Deltadraw[,(1:nvar)],type="l")
abline(h=td[1:nvar],col=(1:nvar))
matplot(out$Deltadraw[,((nvar+1):(2*nvar))],type="l")
abline(h=td[(nvar+1):(2*nvar)],col=(1:nvar))
matplot(out$Vbetadraw[,c(1,7,13,19,25)],type="l")
abline(h=1.5)
matplot(out$Vbetadraw[,-c(1,7,13,19,25)],type="l")
abline(h=.5)
}

}
\seealso{ \code{\link{rhierMnlRwMixture}} }
\keyword{ models}

