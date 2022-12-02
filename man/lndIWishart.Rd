\name{lndIWishart}
\alias{lndIWishart}
\concept{Inverted Wishart distribution}
\concept{density}

\title{Compute Log of Inverted Wishart Density}

\description{
  \code{lndIWishart} computes the log of an Inverted Wishart density.
}

\usage{lndIWishart(nu, V, IW)}

\arguments{
  \item{nu }{ d.f. parameter }
  \item{V  }{ "location" parameter  }
  \item{IW }{ ordinate for density evaluation }
}

\details{
  \eqn{Z} \eqn{\sim}{~} Inverted Wishart(nu,V). \cr
   In this parameterization, \eqn{E[Z] = 1/(nu-k-1) V}, where \eqn{V} is a \eqn{k x k} matrix \cr
  \code{lndIWishart} computes the complete log-density, including normalizing constants.
}

\value{Log density value}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 2, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\seealso{ \code{\link{rwishart}} }

\examples{
lndIWishart(5, diag(3), diag(3)+0.5)
}

\keyword{ distribution }


