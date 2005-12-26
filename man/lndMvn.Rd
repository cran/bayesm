\name{lndMvn}
\alias{lndMvn}
\concept{multivariate normal distribution}
\concept{density}

\title{ Compute Log of Multivariate Normal Density }
\description{
  \code{lndMvn} computes the log of a Multivariate Normal Density.

}

\usage{
lndMvn(x, mu, rooti)
}

\arguments{
  \item{x}{ density ordinate }
  \item{mu}{ mu vector }
  \item{rooti}{ inv of Upper Triangular Cholesky root of Sigma }
}
\details{
  \eqn{z} \eqn{\sim}{~} \eqn{N(mu,\Sigma)}
}
\value{
  log density value
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 2. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{lndMvst}} }
\examples{
##
Sigma=matrix(c(1,.5,.5,1),ncol=2)
lndMvn(x=c(rep(0,2)),mu=c(rep(0,2)),rooti=backsolve(chol(Sigma),diag(2)))
}
\keyword{ distribution }
