\name{rtrun}
\alias{rtrun}
\concept{truncated normal}
\concept{simulation}

\title{ Draw from Truncated Univariate Normal }
\description{
  \code{rtrun} draws from a truncated univariate normal distribution
}
\usage{
rtrun(mu, sigma, a, b)
}
\arguments{
  \item{mu}{ mean }
  \item{sigma}{ sd }
  \item{a}{ lower bound }
  \item{b}{ upper bound }
}
\details{
  Note that due to the vectorization of the rnorm,qnorm commands in R, all arguments can be vectors of 
  equal length.  This makes the inverse CDF method the most efficient to use in R.
}
\value{
  draw (possibly a vector)
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 2. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}
\examples{
##
set.seed(66)
rtrun(mu=c(rep(0,10)),sigma=c(rep(1,10)),a=c(rep(0,10)),b=c(rep(2,10)))
}
\keyword{ distribution }
