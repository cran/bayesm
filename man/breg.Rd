\name{breg}
\alias{breg}
\concept{bayes}
\concept{regression}

\title{Posterior Draws from a Univariate Regression with Unit Error Variance}
\description{
  \code{breg} makes one draw from the posterior of a univariate regression
  (scalar dependent variable) given the error variance = 1.0.
  A natural conjugate, normal prior is used.
}
\usage{
breg(y, X, betabar, A)
}

\arguments{
  \item{y}{ vector of values of dep variable. }
  \item{X}{ n (length(y)) x k Design matrix. }
  \item{betabar}{ k x 1 vector. Prior mean of regression coefficients. }
  \item{A}{ Prior precision matrix. }
}

\details{
  model: \eqn{y=x'\beta + e}.   \eqn{e} \eqn{\sim}{~} \eqn{N(0,1)}.  \cr

  prior: \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})}.
}

\value{
  k x 1 vector containing a draw from the posterior distribution.
}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.

  In particular, X must be a matrix. If you have a vector for X, coerce it
  into a matrix with one column
}

\examples{
##

if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=1000} else {R=10}

## simulate data
set.seed(66)
n=100
X=cbind(rep(1,n),runif(n)); beta=c(1,2)
y=X\%*\%beta+rnorm(n)
##
## set prior
A=diag(c(.05,.05)); betabar=c(0,0)
##
## make draws from posterior
betadraw=matrix(double(R*2),ncol=2)
for (rep in 1:R) {betadraw[rep,]=breg(y,X,betabar,A)}
##
## summarize draws
mat=apply(betadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(beta,mat); rownames(mat)[1]="beta"; print(mat)
}
\keyword{models}
\keyword{regression}
\keyword{distribution}
