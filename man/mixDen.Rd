\name{mixDen}
\alias{mixDen}
\concept{normal mixture}
\concept{marginal distribution}
\concept{density}

\title{ Compute Marginal Density for Multivariate Normal Mixture }
\description{
  \code{mixDen} computes the marginal density for each component of
  a normal mixture at each of the points on a user-specifed grid.
}
\usage{
mixDen(x, pvec, comps)
}
\arguments{
  \item{x}{ array - ith column gives grid points for ith variable }
  \item{pvec}{ vector of mixture component probabilites }
  \item{comps}{ list of lists of components for normal mixture }
}
\details{
  length(comps) is the number of mixture components.  comps[[j]] is a list of
  parameters of the jth component. comps[[j]]$mu is mean vector; comps[[j]]$rooti
  is the UL decomp of \eqn{\Sigma^{-1}}.
}

\value{
  an array of the same dimension as grid with density values.
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 3. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{rnmixGibbs}}  }

\examples{
\dontrun{
##
##  see examples in rnmixGibbs documentation
##
}
}

\keyword{ models }
\keyword{ multivariate }
