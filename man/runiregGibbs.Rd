\name{runiregGibbs}
\alias{runiregGibbs}
\concept{bayes}
\concept{Gibbs Sampler}
\concept{regression}
\concept{MCMC}

\title{ Gibbs Sampler for Univariate Regression }
\description{
 \code{runiregGibbs} implements a Gibbs Sampler to draw from posterior of a univariate regression with a conditionally conjugate prior.
}
\usage{
runiregGibbs(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y,X)}
  \item{Prior}{ list(betabar,A, nu, ssq) }
  \item{Mcmc}{ list(sigmasq,R,keep)}
}
\details{
  Model: \eqn{y = Xbeta + e}.  \eqn{e} \eqn{\sim}{~} \eqn{N(0,sigmasq)}. \cr

  Priors: \eqn{beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})}.  \eqn{sigmasq} \eqn{\sim}{~} \eqn{(nu*ssq)/chisq_{nu}}.
  List arguments contain  
  \itemize{
    \item{\code{X}}{n x k Design Matrix}
    \item{\code{y}}{n x 1 vector of observations} 
    \item{\code{betabar}}{k x 1 prior mean (def: 0)}
    \item{\code{A}}{k x k prior precision matrix (def: .01I)} 
    \item{\code{nu}}{ d.f. parm for Inverted Chi-square prior (def: 3)}
    \item{\code{ssq}}{ scale parm for Inverted Chi-square prior (def:var(y))}
    \item{\code{R}}{ number of MCMC draws }
    \item{\code{keep}}{ thinning parameter - keep every keepth draw }
  }
}
\value{
  list of MCMC draws
  \item{betadraw }{ R x k array of betadraws }
  \item{sigmasqdraw }{ R vector of sigma-sq draws}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 3. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\seealso{ \code{\link{runireg}} }
\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=1000} else {R=10}
set.seed(66)
n=100
X=cbind(rep(1,n),runif(n)); beta=c(1,2); sigsq=.25
y=X\%*\%beta+rnorm(n,sd=sqrt(sigsq))

A=diag(c(.05,.05)); betabar=c(0,0)
nu=3; ssq=1.0

Data=list(y=y,X=X); Mcmc=list(R=R,keep=1) ; Prior=list(A=A,betabar=betabar,nu=nu,ssq=ssq)
out=runiregGibbs(Data=Data,Prior=Prior,Mcmc=Mcmc)
cat(" betadraws ",fill=TRUE)
mat=apply(out$betadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(beta,mat); rownames(mat)[1]="beta"; print(mat)
cat(" Sigma-sq draws",fill=TRUE)
cat(" sigma-sq= ",sigsq,fill=TRUE)
print(quantile(out$sigmasqdraw,probs=c(.01,.05,.5,.95,.99)))
}
\keyword{ regression }
