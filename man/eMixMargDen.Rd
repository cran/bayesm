\name{eMixMargDen}
\alias{eMixMargDen}
\concept{normal mixtures}
\concept{bayes}
\concept{MCMC}

\title{ Compute Marginal Densities of A Normal Mixture Averaged over MCMC Draws }
\description{
  \code{eMixMargDen} assumes that a multivariate mixture of normals has been fitted
  via MCMC (using \code{rnmixGibbs}).  For each MCMC draw, the marginal densities
  for each component in the multivariate mixture are computed on a user-supplied
  grid and then averaged over draws. 
}

\usage{
eMixMargDen(grid, probdraw, compdraw)
}

\arguments{
  \item{grid}{  array of grid points, grid[,i] are ordinates for ith dimension of the density }
  \item{probdraw}{ array - each row of which contains a draw of probabilities of mixture comp }
  \item{compdraw}{ list of lists of draws of mixture comp moments }
}

\details{
  length(compdraw) is number of MCMC draws. \cr
  compdraw[[i]] is a list draws of mu and inv Chol root for each of mixture components. \cr
  compdraw[[i]][[j]] is jth component.  compdraw[[i]][[j]]$mu is mean vector; compdraw[[i]][[j]]$rooti
  is the UL decomp of \eqn{\Sigma^{-1}}.
}

\value{
  an array of the same dimension as grid with density values.
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type. To avoid errors, call with
  output from \code{\link{rnmixGibbs}}.
}

\seealso{ \code{\link{rnmixGibbs}}  }

\keyword{ models }
\keyword{ multivariate }
