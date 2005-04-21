\name{lndIChisq}
\alias{lndIChisq}
\concept{Inverted Chi-squared Distribution}
\concept{density}

\title{ Compute Log of Inverted Chi-Squared Density }
\description{
  \code{lndIChisq} computes the log of an Inverted Chi-Squared Density.
}
\usage{
lndIChisq(nu, ssq, x)
}
\arguments{
  \item{nu}{ d.f. parameter }
  \item{ssq}{ scale parameter }
  \item{x}{ ordinate for density evaluation }
}
\details{
  \eqn{Z= \nu*ssq/\chi^2_{\nu}}, \eqn{Z} \eqn{\sim}{~} Inverted Chi-Squared.  \cr
  \code{lndIChisq} computes the complete log-density, including normalizing constants.
}
\value{
  log density value
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 2. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{dchisq}} }
\examples{
##
lndIChisq(3,1,2)
}
\keyword{ distribution }
