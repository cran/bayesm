\name{mixDenBi}
\alias{mixDenBi}
\concept{normal mixture}
\concept{marginal distribution}
\concept{density}

\title{Compute Bivariate Marginal Density for a Normal Mixture}

\description{
\code{mixDenBi} computes the implied bivariate marginal density from a mixture of normals with specified mixture probabilities and component parameters. 
}

\usage{mixDenBi(i, j, xi, xj, pvec, comps)}

\arguments{
  \item{i}{ index of first variable }
  \item{j}{ index of second variable }
  \item{xi}{ grid of values of first variable }
  \item{xj}{ grid of values of second variable }
  \item{pvec}{ normal mixture probabilities }
  \item{comps}{ list of lists of components }
}

\details{
  \tabular{ll}{
    \code{length(comps)   } \tab is the number of mixture components \cr
    \code{comps[[j]]      } \tab is a list of parameters of the \eqn{j}th component \cr
    \code{comps[[j]]$mu   } \tab is mean vector \cr
    \code{comps[[j]]$rooti} \tab is the UL decomp of \eqn{\Sigma^{-1}}
  }
}

\value{An array (\code{length(xi)=length(xj) x 2}) with density value}
  
\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{ Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\seealso{ \code{\link{rnmixGibbs}}, \code{\link{mixDen}}  }

\keyword{models}
\keyword{multivariate}

