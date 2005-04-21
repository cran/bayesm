\name{numEff}
\alias{numEff}
\concept{numerical efficiency}

\title{ Compute Numerical Standard Error and Relative Numerical Efficiency }
\description{
  \code{numEff} computes the numerical standard error for the mean of a vector of draws as well as the relative
  numerical efficiency (ratio of variance of mean of this time series process relative to iid sequence).
}

\usage{
numEff(x, m = max(length(x), (100/sqrt(5000)) * sqrt(length(x))))
}

\arguments{
  \item{x}{ R x 1 vector of draws }
  \item{m}{ number of lags for autocorrelations }
}

\details{
  default for number of lags is chosen so that if R = 5000, m =100 and increases as the sqrt(R).
}
\value{
  \item{stderr }{standard error of the mean of x}
  \item{f }{ variance ratio (relative numerical efficiency) }
}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 3. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

}
\keyword{ ts }
\keyword{ utilities }
