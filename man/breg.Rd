\name{breg}
\alias{breg}

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
  by Allenby, McCulloch, and Rossi. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\examples{
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
R=1000
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
