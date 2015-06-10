\name{runireg}
\alias{runireg}
\concept{bayes}
\concept{regression}

\title{ IID Sampler for Univariate Regression }
\description{
 \code{runireg} implements an iid sampler to draw from the posterior of a univariate regression with a conjugate prior.
}
\usage{
runireg(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y,X)}
  \item{Prior}{ list(betabar,A, nu, ssq) }
  \item{Mcmc}{ list(R,keep,nprint)}
}
\details{
  Model: \eqn{y = X\beta + e}.  \eqn{e} \eqn{\sim}{~} \eqn{N(0,\sigma^2)}. \cr

  Priors: \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar,\sigma^2*A^{-1})}. 
 \eqn{\sigma^2} \eqn{\sim}{~} \eqn{(nu*ssq)/\chi^2_{nu}}.
  List arguments contain  
  \itemize{
    \item{\code{X}}{n x k Design Matrix}
    \item{\code{y}}{n x 1 vector of observations} 
    \item{\code{betabar}}{k x 1 prior mean (def: 0)}
    \item{\code{A}}{k x k prior precision matrix (def: .01I)} 
    \item{\code{nu}}{ d.f. parm for Inverted Chi-square prior (def: 3)}
    \item{\code{ssq}}{ scale parm for Inverted Chi-square prior (def: var(y))}
    \item{\code{R}}{ number of draws }
    \item{\code{keep}}{ thinning parameter - keep every keepth draw}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
  }
}
\value{
  list of iid draws
  \item{betadraw }{ R x k array of betadraws }
  \item{sigmasqdraw }{ R vector of sigma-sq draws}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 2. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\seealso{ \code{\link{runiregGibbs}} }
\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)
n=200
X=cbind(rep(1,n),runif(n)); beta=c(1,2); sigsq=.25
y=X\%*\%beta+rnorm(n,sd=sqrt(sigsq))

out=runireg(Data=list(y=y,X=X),Mcmc=list(R=R))

cat("Summary of beta/sigma-sq draws",fill=TRUE)
summary(out$betadraw,tvalues=beta)
summary(out$sigmasqdraw,tvalues=sigsq)

if(0){
## plotting examples
plot(out$betadraw)
}

}
\keyword{ regression }
