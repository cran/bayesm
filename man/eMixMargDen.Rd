\name{eMixMargDen}
\alias{eMixMargDen}
\concept{normal mixtures}
\concept{bayes}
\concept{MCMC}

\title{Compute Marginal Densities of A Normal Mixture Averaged over MCMC Draws}

\description{
  \code{eMixMargDen} assumes that a multivariate mixture of normals has been fitted
  via MCMC (using \code{rnmixGibbs}). For each MCMC draw, \code{eMixMargDen} computes 
  the marginal densities for each component in the multivariate mixture on a user-supplied
  grid and then averages over the MCMC draws. 
}

\usage{eMixMargDen(grid, probdraw, compdraw)}

\arguments{
  \item{grid}{ array of grid points, \code{grid[,i]} are ordinates for \code{i}th dimension of the density }
  \item{probdraw}{ array where each row contains a draw of probabilities of the mixture component }
  \item{compdraw}{ list of lists of draws of mixture component moments }
}

\details{
  \tabular{ll}{
    \code{length(compdraw)        }  \tab is the number of MCMC draws. \cr
    \code{compdraw[[i]]           }  \tab is a list draws of mu and of the inverse Cholesky root for each of mixture components. \cr
    \code{compdraw[[i]][[j]]      }  \tab is \code{j}th component. \cr
    \code{compdraw[[i]][[j]]$mu   }  \tab is mean vector. \cr
    \code{compdraw[[i]][[j]]$rooti} \tab is the UL decomp of \eqn{\Sigma^{-1}}.
  }
}

\value{An array of the same dimension as \code{grid} with density values.}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type. To avoid errors, call with output from \code{\link{rnmixGibbs}}.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}by Rossi, Allenby, and McCulloch. }

\seealso{\code{\link{rnmixGibbs}}}

\keyword{models}
\keyword{multivariate}
