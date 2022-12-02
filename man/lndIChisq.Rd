\name{lndIChisq}
\alias{lndIChisq}
\concept{Inverted Chi-squared Distribution}
\concept{density}

\title{Compute Log of Inverted Chi-Squared Density}

\description{
\code{lndIChisq} computes the log of an Inverted Chi-Squared Density.
}

\usage{lndIChisq(nu, ssq, X)}

\arguments{
  \item{nu  }{ d.f. parameter }
  \item{ssq }{ scale parameter }
  \item{X   }{ ordinate for density evaluation (this must be a matrix)}
}

\details{
  \eqn{Z = nu*ssq / \chi^2_{nu}} with \eqn{Z} \eqn{\sim}{~} Inverted Chi-Squared.  \cr
  \code{lndIChisq} computes the complete log-density, including normalizing constants.
}

\value{Log density value}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 2, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{dchisq}} }

\examples{
lndIChisq(3, 1, matrix(2))
}

\keyword{distribution}
