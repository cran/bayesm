\name{lndIWishart}
\alias{lndIWishart}
\concept{Inverted Wishart distribution}
\concept{density}

\title{ Compute Log of Inverted Wishart Density }
\description{
  \code{lndIWishart} computes the log of an Inverted Wishart density.
}
\usage{
lndIWishart(nu, S, IW)
}
\arguments{
  \item{nu}{ d.f. parameter }
  \item{S}{ "location" parameter  }
  \item{IW}{ ordinate for density evaluation }
}
\details{
  \eqn{Z = Wishart(nu,V^{-1})^{-1}}, \eqn{Z} \eqn{\sim}{~} Inverted Wishart(nu,V).    \cr
  \code{lndIWishart} computes the complete log-density, including normalizing constants.
}
\value{
  log density value
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 2. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{rwishart}} }
\examples{
##
lndIWishart(5,diag(3),(diag(3)+.5))
}
\keyword{ distribution }


