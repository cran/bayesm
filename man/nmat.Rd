\name{nmat}
\alias{nmat}

\title{Convert Covariance Matrix to a Correlation Matrix}

\description{
\code{nmat} converts a covariance matrix (stored as a vector, col by col) to a correlation matrix (also stored as a vector).
}

\usage{nmat(vec)}

\arguments{
  \item{vec}{ \eqn{k x k} Cov matrix stored as a \eqn{k*k x 1} vector (col by col) }
}

\details{
This routine is often used with apply to convert an \eqn{R x (k*k)} array of covariance MCMC draws to correlations. As in \code{corrdraws = apply(vardraws, 1, nmat)}.
}

\value{\eqn{k*k x 1} vector with correlation matrix}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\examples{
set.seed(66)
X = matrix(rnorm(200,4), ncol=2)
Varmat = var(X)
nmat(as.vector(Varmat))
}

\keyword{utilities}
\keyword{array}

