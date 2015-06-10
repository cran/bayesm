\name{rmixGibbs}
\alias{rmixGibbs}
\title{ Gibbs Sampler for Normal Mixtures w/o Error Checking}
\description{
  \code{rmixGibbs} makes one draw using the Gibbs Sampler for a mixture of multivariate normals.
}
\usage{
rmixGibbs(y, Bbar, A, nu, V, a, p, z)
}
\arguments{
  \item{y}{ data array - rows are obs }
  \item{Bbar}{ prior mean for mean vector of each norm comp }
  \item{A}{ prior precision parameter}
  \item{nu}{ prior d.f. parm }
  \item{V}{ prior location matrix for covariance priro }
  \item{a}{ Dirichlet prior parms }
  \item{p}{ prior prob of each mixture component }
  \item{z}{ component identities for each observation -- "indicators"}
}
\details{
  \code{rmixGibbs} is not designed to be called directly. Instead, use \code{rnmixGibbs} wrapper function.
}
\value{
  a list containing:
  \item{p}{draw mixture probabilities }
  \item{z}{draw of indicators of each component}
  \item{comps}{new draw of normal component parameters }
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 5. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Rob McCulloch and Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{rnmixGibbs}} }

\keyword{ multivariate } 

