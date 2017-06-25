\name{rmixture}
\alias{rmixture}
\concept{mixture of normals}
\concept{simulation}

\title{Draw from Mixture of Normals}

\description{
  \code{rmixture} simulates iid draws from a Multivariate Mixture of Normals
}

\usage{rmixture(n, pvec, comps)}

\arguments{
  \item{n     }{ number of observations }
  \item{pvec  }{ \eqn{ncomp x 1} vector of prior probabilities for each mixture component }
  \item{comps }{ list of mixture component parameters }
}

\details{
  \code{comps} is a list of length \code{ncomp} with \code{ncomp = length(pvec)}. \cr
  \code{comps[[j]][[1]]} is mean vector for the \eqn{j}th component. \cr
  \code{comps[[j]][[2]]} is the inverse of the cholesky root of \eqn{\Sigma} for \eqn{j}th component
}

\value{
  A list containing:
  \item{x: }{ an \eqn{n x} \code{length(comps[[1]][[1]])} array of iid draws }
  \item{z: }{ an \eqn{n x 1} vector of indicators of which component each draw is taken from }
}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\seealso{ \code{\link{rnmixGibbs}} }

\keyword{multivariate}
\keyword{distribution}
