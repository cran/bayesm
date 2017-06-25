\name{rtrun}
\alias{rtrun}
\concept{truncated normal}
\concept{simulation}

\title{Draw from Truncated Univariate Normal}

\description{
\code{rtrun} draws from a truncated univariate normal distribution.
}

\usage{rtrun(mu, sigma, a, b)}

\arguments{
  \item{mu}{ mean }
  \item{sigma}{ standard deviation }
  \item{a}{ lower bound }
  \item{b}{ upper bound }
}

\details{
Note that due to the vectorization of the \code{rnorm} and \code{qnorm} commands in R, all arguments can be vectors of equal length. This makes the inverse CDF method the most efficient to use in R.
}

\value{Draw (possibly a vector)}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.

**Note also that \code{rtrun} is currently affected by the numerical accuracy of the inverse CDF function when trunctation points are far out in the tails of the distribution, where \dQuote{far out} means \eqn{|a - \mu| / \sigma > 6} and/or \eqn{|b - \mu| / \sigma > 6}. A fix will be implemented in a future version of \code{bayesm}.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 2, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. \cr \url{http://www.perossi.org/home/bsm-1}}

\examples{
set.seed(66)
rtrun(mu=c(rep(0,10)), sigma=c(rep(1,10)), a=c(rep(0,10)), b=c(rep(2,10)))
}

\keyword{distribution}
