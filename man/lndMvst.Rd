\name{lndMvst}
\alias{lndMvst}
\concept{multivariate t distribution}
\concept{student-t distribution}
\concept{density}

\title{Compute Log of Multivariate Student-t Density}

\description{
  \code{lndMvst} computes the log of a Multivariate Student-t Density.
}

\usage{lndMvst(x, nu, mu, rooti, NORMC)}

\arguments{
  \item{x    }{ density ordinate }
  \item{nu   }{ d.f. parameter }
  \item{mu   }{ mu vector }
  \item{rooti}{ inv of Cholesky root of \eqn{\Sigma} }
  \item{NORMC}{ include normalizing constant (def: \code{FALSE}) }
}

\details{
  \eqn{z} \eqn{\sim}{~} \eqn{MVst(mu, nu, \Sigma)}
}

\value{Log density value}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 2, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\seealso{ \code{\link{lndMvn}} }

\examples{
Sigma = matrix(c(1, 0.5, 0.5, 1), ncol=2)
lndMvst(x=c(rep(0,2)), nu=4,mu=c(rep(0,2)), rooti=backsolve(chol(Sigma),diag(2)))
}

\keyword{distribution}
