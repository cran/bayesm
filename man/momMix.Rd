\name{momMix}
\alias{momMix}
\concept{mcmc}
\concept{normal mixture}
\concept{posterior moments}

\title{ Compute Posterior Expectation of Normal Mixture Model Moments }
\description{
  \code{momMix} averages the moments of a normal mixture model over MCMC draws.
}
\usage{
momMix(probdraw, compdraw)
}

\arguments{
  \item{probdraw}{ R x ncomp list of draws of mixture probs }
  \item{compdraw}{ list of length R of draws of mixture component moments }
}
\details{
  R is the number of MCMC draws in argument list above. \cr
  ncomp is the number of mixture components fitted.\cr
  compdraw is a list of lists of lists with mixture components.  \cr
  compdraw[[i]] is ith draw. \cr
  compdraw[[i]][[j]][[1]] is the mean parameter vector for the jth component, ith MCMC draw. \cr
  compdraw[[i]][[j]][[2]] is the UL decomposition of \eqn{\Sigma^{-1}} for the jth component, ith MCMC draw. 
 
}
\value{
  a list of the following items \dots
  \item{mu }{Posterior Expectation of Mean}
  \item{sigma }{Posterior Expecation of Covariance Matrix}
  \item{sd }{Posterior Expectation of Vector of Standard Deviations}
  \item{corr }{Posterior Expectation of Correlation Matrix}
}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 5. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}
\seealso{ \code{\link{rmixGibbs}}}

\keyword{ multivariate }
