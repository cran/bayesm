\name{lndMvst}
\alias{lndMvst}

\title{ Compute Log of Multivariate Student-t Density }
\description{
  \code{lndIChisq} computes the log of a Multivariate Student-t Density.
}
\usage{
lndMvst(x, nu, mu, rooti)
}

\arguments{
  \item{x}{ density ordinate }
  \item{nu}{ d.f. parameter }
  \item{mu}{ mu vector }
  \item{rooti}{ inv of Cholesky root of Sigma }
}

\details{
  \eqn{z} \eqn{\sim}{~} \eqn{MVst(mu,nu,\Sigma)}
  note: does not include full normalizing constant 
}
\value{
  log density value
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 2. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{lndMvn}} }
\examples{
##
Sigma=matrix(c(1,.5,.5,1),ncol=2)
lndMvst(x=c(rep(0,2)),nu=4,mu=c(rep(0,2)),rooti=backsolve(chol(Sigma),diag(2)))
}
\keyword{ distribution }
