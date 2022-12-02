\name{condMom}
\alias{condMom}
\concept{normal distribution}
\concept{conditional distribution}

\title{Computes Conditional Mean/Var of One Element of MVN given All Others}

\description{
\code{condMom} compute moments of conditional distribution of the \eqn{i}th element of a multivariate normal given all others.
}

\usage{condMom(x, mu, sigi, i)}

\arguments{
  \item{x}{ vector of values to condition on; \eqn{i}th element not used }
  \item{mu}{ mean vector with \code{length(x)} = \eqn{n} }
  \item{sigi}{ inverse of covariance matrix; dimension \eqn{n x n} }
  \item{i}{ conditional distribution of \eqn{i}th element }
}

\details{
  \eqn{x} \eqn{\sim}{~} \eqn{MVN(mu, sigi^{-1})}.

  \code{condMom} computes moments of \eqn{x_i} given \eqn{x_{-i}}.
}

\value{
  A list containing:
  \item{cmean }{ conditional mean }
  \item{cvar  }{ conditional variance}
}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\examples{
sig  = matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), ncol=3)
sigi = chol2inv(chol(sig))
mu   = c(1,2,3)
x    = c(1,1,1)

condMom(x, mu, sigi, 2)
}

\keyword{distribution}
