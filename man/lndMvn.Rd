\name{lndMvn}
\alias{lndMvn}
\concept{multivariate normal distribution}
\concept{density}

\title{ Compute Log of Multivariate Normal Density }

\description{
\code{lndMvn} computes the log of a Multivariate Normal Density.
}

\usage{lndMvn(x, mu, rooti)}

\arguments{
  \item{x    }{ density ordinate }
  \item{mu   }{ mu vector }
  \item{rooti}{ inv of upper triangular Cholesky root of \eqn{\Sigma} }
}

\details{
  \eqn{z} \eqn{\sim}{~} \eqn{N(mu,\Sigma)}
}

\value{Log density value}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 2, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\seealso{ \code{\link{lndMvst}} }

\examples{
Sigma = matrix(c(1, 0.5, 0.5, 1), ncol=2)
lndMvn(x=c(rep(0,2)), mu=c(rep(0,2)), rooti=backsolve(chol(Sigma),diag(2)))
}

\keyword{distribution}
