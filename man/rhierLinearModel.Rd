\name{rhierLinearModel}
\alias{rhierLinearModel}
\concept{bayes}
\concept{MCMC}
\concept{Gibbs Sampling}
\concept{hierarchical models}
\concept{linear model}

\title{ Gibbs Sampler for Hierarchical Linear Model }
\description{
  \code{rhierLinearModel} implements a Gibbs Sampler for hierarchical linear models.
}
\usage{
rhierLinearModel(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(regdata,Z) (Z optional). }
  \item{Prior}{ list(Deltabar,A,nu.e,ssq,nu,V)  (optional).}
  \item{Mcmc}{ list(R,keep) (R required).}
}
\details{
  Model: length(regdata) regression equations. \cr
        \eqn{y_i = X_ibeta_i + e_i}. \eqn{e_i} \eqn{\sim}{~} \eqn{N(0,tau_i)}.  nvar X vars in each equation. 

 Priors:\cr
        \eqn{tau_i} \eqn{\sim}{~} nu.e*\eqn{ssq_i/\chi^2_{nu.e}}.  \eqn{tau_i} is the variance of \eqn{e_i}.\cr
        \eqn{beta_i} \eqn{\sim}{~} N(ZDelta[i,],\eqn{V_{beta}}). \cr
               Note:  ZDelta is the matrix Z * Delta; [i,] refers to ith row of this product.

          \eqn{vec(Delta)} given \eqn{V_{beta}} \eqn{\sim}{~} \eqn{N(vec(Deltabar),V_{beta} (x) A^{-1})}.\cr
          \eqn{V_{beta}} \eqn{\sim}{~} \eqn{IW(nu,V)}. \cr
              \eqn{Delta, Deltabar} are nz x nvar.  \eqn{A} is nz x nz.  \eqn{V_{beta}} is nvar x nvar.
        
          Note: if you don't have any z vars, set Z=iota (nreg x 1).

  List arguments contain:
  \itemize{
    \item{\code{regdata}}{ list of lists with X,y matrices for each of length(regdata) regressions}
    \item{\code{regdata[[i]]$X}}{ X matrix for equation i }
    \item{\code{regdata[[i]]$y}}{ y vector for equation i }
    \item{\code{Deltabar}}{ nz x nvar matrix of prior means (def: 0)}
    \item{\code{A}}{ nz x nz matrix for prior precision (def: .01I)}
    \item{\code{nu.e}}{ d.f. parm for regression error variance prior (def: 3)}
    \item{\code{ssq}}{ scale parm for regression error var prior (def: var(\eqn{y_i}))}
    \item{\code{nu}}{ d.f. parm for Vbeta prior (def: nvar+3)}
    \item{\code{V}}{ Scale location matrix for Vbeta prior (def: nu*I)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
   }
}
\value{
  a list containing
  \item{betadraw}{nreg x nvar x R/keep array of individual regression coef draws}
  \item{taudraw}{R/keep x nreg array of error variance draws}
  \item{Deltadraw}{R/keep x nz x nvar array of Deltadraws}
  \item{Vbetadraw}{R/keep x nvar*nvar array of Vbeta draws}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 3. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) # set env var LONG_TEST to run
{

nreg=100; nobs=100; nvar=3
Vbeta=matrix(c(1,.5,0,.5,2,.7,0,.7,1),ncol=3)
Z=cbind(c(rep(1,nreg)),3*runif(nreg)); Z[,2]=Z[,2]-mean(Z[,2])
nz=ncol(Z)
Delta=matrix(c(1,-1,2,0,1,0),ncol=2)
Delta=t(Delta) # first row of Delta is means of betas
Beta=matrix(rnorm(nreg*nvar),nrow=nreg)\%*\%chol(Vbeta)+Z\%*\%Delta
tau=.1
iota=c(rep(1,nobs))
regdata=NULL
for (reg in 1:nreg) { X=cbind(iota,matrix(runif(nobs*(nvar-1)),ncol=(nvar-1)))
	y=X\%*\%Beta[reg,]+sqrt(tau)*rnorm(nobs); regdata[[reg]]=list(y=y,X=X) }

nu.e=3
ssq=NULL
for(reg in 1:nreg) {ssq[reg]=var(regdata[[reg]]$y)}
nu=nvar+3
V=nu*diag(c(rep(1,nvar)))
A=diag(c(rep(.01,nz)),ncol=nz)
Deltabar=matrix(c(rep(0,nz*nvar)),nrow=nz)

R=2000

Data=list(regdata=regdata,Z=Z)
Prior=list(Deltabar=Deltabar,A=A,nu.e=nu.e,ssq=ssq,nu=nu,V=V)
Mcmc=list(R=R,keep=1)
out=rhierLinearModel(Data=Data,Mcmc=Mcmc)

cat(" Deltadraws ",fill=TRUE)
mat=apply(out$Deltadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(as.vector(Delta),mat); rownames(mat)[1]="delta"; print(mat)
cat(" Vbetadraws ",fill=TRUE)
mat=apply(out$Vbetadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(as.vector(Vbeta),mat); rownames(mat)[1]="Vbeta"; print(mat)

}

}
\keyword{ regression }
