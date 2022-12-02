\name{momMix}
\alias{momMix}
\concept{mcmc}
\concept{normal mixture}
\concept{posterior moments}

\title{Compute Posterior Expectation of Normal Mixture Model Moments}

\description{
\code{momMix} averages the moments of a normal mixture model over MCMC draws.
}

\usage{momMix(probdraw, compdraw)}

\arguments{
  \item{probdraw}{ \eqn{R x ncomp} list of draws of mixture probs }
  \item{compdraw}{ list of length \eqn{R} of draws of mixture component moments }
}

\details{
  \tabular{ll}{
    \code{R                       } \tab is the number of MCMC draws in argument list above. \cr
    \code{ncomp                   } \tab is the number of mixture components fitted.\cr
    \code{compdraw                } \tab is a list of lists of lists with mixture components.  \cr
    \code{compdraw[[i]]           } \tab is \eqn{i}th draw. \cr
    \code{compdraw[[i]][[j]][[1]] } \tab is the mean parameter vector for the \eqn{j}th component, \eqn{i}th MCMC draw. \cr
    \code{compdraw[[i]][[j]][[2]] } \tab is the UL decomposition of \eqn{\Sigma^{-1}} for the \eqn{j}th component, \eqn{i}th MCMC draw
  }
}
\value{
  A list containing:
  \item{mu    }{posterior expectation of mean}
  \item{sigma }{posterior expectation of covariance matrix}
  \item{sd    }{posterior expectation of vector of standard deviations}
  \item{corr  }{posterior expectation of correlation matrix}
}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\seealso{ \code{\link{rmixGibbs}}}

\keyword{multivariate}
