\name{mixDen}
\alias{mixDen}
\concept{normal mixture}
\concept{marginal distribution}
\concept{density}

\title{Compute Marginal Density for Multivariate Normal Mixture}

\description{
\code{mixDen} computes the marginal density for each dimension of a normal mixture at each of the points on a user-specifed grid.}

\usage{mixDen(x, pvec, comps)}

\arguments{
  \item{x}{ array where \eqn{i}th column gives grid points for \eqn{i}th variable }
  \item{pvec}{ vector of mixture component probabilites }
  \item{comps}{ list of lists of components for normal mixture }
}

\details{
  \tabular{ll}{
    \code{length(comps)   } \tab is the number of mixture components \cr
    \code{comps[[j]]      } \tab is a list of parameters of the \eqn{j}th component \cr
    \code{comps[[j]]$mu   } \tab is mean vector \cr
    \code{comps[[j]]$rooti} \tab is the UL decomp of \eqn{\Sigma^{-1}}
  }
}

\value{An array of the same dimension as grid with density values}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\seealso{ \code{\link{rnmixGibbs}} }

\keyword{models}
\keyword{multivariate}
