\name{rmixGibbs}
\alias{rmixGibbs}
\title{ Gibbs Sampler for Normal Mixtures w/o Error Checking}
\description{
  \code{rmixGibbs} makes one draw using the Gibbs Sampler for a mixture of multivariate normals. \code{rmixGibbs} is not designed to be called directly. Instead, use \code{rnmixGibbs} wrapper function.
}
\usage{
rmixGibbs(y, Bbar, A, nu, V, a, p, z)
}
\arguments{
  \item{y    }{ data array where rows are obs }
  \item{Bbar }{ prior mean for mean vector of each norm comp }
  \item{A    }{ prior precision parameter}
  \item{nu   }{ prior d.f. parm }
  \item{V    }{ prior location matrix for covariance prior }
  \item{a    }{ Dirichlet prior parms }
  \item{p    }{ prior prob of each mixture component }
  \item{z    }{ component identities for each observation -- "indicators" }
}

\value{
  a list containing:
  \item{p}{ draw of mixture probabilities}
  \item{z}{ draw of indicators of each component}
  \item{comps}{ new draw of normal component parameters }
}

\author{Rob McCulloch (Arizona State University) and Peter Rossi (Anderson School, UCLA), \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 5 \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\seealso{ \code{\link{rnmixGibbs}} }

\keyword{multivariate}

