\name{rmixture}
\alias{rmixture}
\concept{mixture of normals}
\concept{simulation}

\title{ Draw from Mixture of Normals }
\description{
  \code{rmixture} simulates iid draws from a Multivariate Mixture of Normals
}
\usage{
rmixture(n, pvec, comps)
}
\arguments{
  \item{n}{ number of observations }
  \item{pvec}{ ncomp x 1 vector of prior probabilities for each mixture component }
  \item{comps}{ list of mixture component parameters }
}
\details{
  comps is a list of length, ncomp = length(pvec). comps[[j]][[1]] is mean vector for the jth component. 
  comps[[j]][[2]] is the inverse of the cholesky root of \eqn{\Sigma} for that component
}
\value{
  A list containing \ldots
  \item{x}{ An n x length(comps[[1]][[1]]) array of iid draws }
  \item{z}{ A n x 1 vector of indicators of which component each draw is taken from }
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{rnmixGibbs}} }
\keyword{ distribution }
\keyword{ multivariate }
