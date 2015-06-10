\name{logMargDenNR}
\alias{logMargDenNR}
\concept{Newton-Raftery approximation}
\concept{bayes}
\concept{marginal likelihood}
\concept{density}

\title{ Compute Log Marginal Density Using Newton-Raftery Approx }
\description{
  \code{logMargDenNR} computes log marginal density using the Newton-Raftery approximation.\cr
  Note: this approximation can be influenced by outliers in the vector of log-likelihoods. Use 
  with \strong{care} . 
}
\usage{
logMargDenNR(ll)
}
\arguments{
  \item{ll}{ vector of log-likelihoods evaluated at length(ll) MCMC draws }
}
\value{
  approximation to log marginal density value.
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 6. \cr
  \url{http://www.perossi.org/home/bsm-1}
}
\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\keyword{ distribution }
