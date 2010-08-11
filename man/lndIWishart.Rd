\name{lndIWishart}
\alias{lndIWishart}
\concept{Inverted Wishart distribution}
\concept{density}

\title{ Compute Log of Inverted Wishart Density }
\description{
  \code{lndIWishart} computes the log of an Inverted Wishart density.
}
\usage{
lndIWishart(nu, V, IW)
}
\arguments{
  \item{nu}{ d.f. parameter }
  \item{V}{ "location" parameter  }
  \item{IW}{ ordinate for density evaluation }
}
\details{
  \eqn{Z} \eqn{\sim}{~} Inverted Wishart(nu,V).    \cr
   in this parameterization, \eqn{E[Z]=1/(nu-k-1) V}, V is a k x k matrix
  \code{lndIWishart} computes the complete log-density, including normalizing constants.
}
\value{
  log density value
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

\seealso{ \code{\link{rwishart}} }
\examples{
##
lndIWishart(5,diag(3),(diag(3)+.5))
}
\keyword{ distribution }


