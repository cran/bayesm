\name{logMargDenNR}
\alias{logMargDenNR}
\concept{Newton-Raftery approximation}
\concept{bayes}
\concept{marginal likelihood}
\concept{density}

\title{Compute Log Marginal Density Using Newton-Raftery Approx}

\description{
\code{logMargDenNR} computes log marginal density using the Newton-Raftery approximation. 
}

\usage{logMargDenNR(ll)}

\arguments{
  \item{ll}{ vector of log-likelihoods evaluated at \code{length(ll)} MCMC draws }
}

\value{Approximation to log marginal density value.}

\section{Warning}{
This approximation can be influenced by outliers in the vector of log-likelihoods; use with \strong{care}.  This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{ Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 6, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }


\keyword{distribution}
