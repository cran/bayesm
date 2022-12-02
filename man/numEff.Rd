\name{numEff}
\alias{numEff}
\concept{numerical efficiency}

\title{Compute Numerical Standard Error and Relative Numerical Efficiency}

\description{
\code{numEff} computes the numerical standard error for the mean of a vector of draws as well as the relative numerical efficiency (ratio of variance of mean of this time series process relative to iid sequence).
}

\usage{numEff(x, m = as.integer(min(length(x),(100/sqrt(5000))*sqrt(length(x)))))}

\arguments{
  \item{x}{ \eqn{R x 1} vector of draws }
  \item{m}{ number of lags for autocorrelations }
}

\details{
default for number of lags is chosen so that if \eqn{R=5000}, \eqn{m=100} and increases as the \eqn{sqrt(R)}.
}

\value{
  A list containing:
  \item{stderr }{standard error of the mean of \eqn{x}}
  \item{f }{ variance ratio (relative numerical efficiency) }
}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
numEff(rnorm(1000), m=20)
numEff(rnorm(1000))
}

\keyword{ts}
\keyword{utilities}
