\name{mixDenBi}
\alias{mixDenBi}
\concept{normal mixture}
\concept{marginal distribution}
\concept{density}

\title{ Compute Bivariate Marginal Density for a Normal Mixture }
\description{
  \code{mixDenBi} computes the implied bivariate marginal density from a mixture of
  normals with specified mixture probabilities and component parameters. 
}
\usage{
mixDenBi(i, j, xi, xj, pvec, comps)
}
\arguments{
  \item{i}{ index of first variable }
  \item{j}{ index of second variable }
  \item{xi}{ grid of values of first variable }
  \item{xj}{ grid of values of second variable }
  \item{pvec}{ normal mixture probabilities }
  \item{comps}{ list of lists of components }
}
\details{
  length(comps) is the number of mixture components.  comps[[j]] is a list of
  parameters of the jth component. comps[[j]]$mu is mean vector; comps[[j]]$rooti
  is the UL decomp of \eqn{\Sigma^{-1}}.
}

\value{
  an array (length(xi)=length(xj) x 2) with density value
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

\seealso{ \code{\link{rnmixGibbs}}, \code{\link{mixDen}}  }

\examples{
\dontrun{
##
##  see examples in rnmixGibbs documentation
##
}
}

\keyword{ models }
\keyword{ multivariate }

