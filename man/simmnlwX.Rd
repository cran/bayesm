\name{simmnlwX}
\alias{simmnlwX}
\title{ Simulate from MNL given X Matrix }
\description{
  \code{simmnlwX} simulates from MNL given X Matrix. 
}
\usage{
simmnlwX(n, X, beta)
}
\arguments{
  \item{n}{ number of obs }
  \item{X}{ n*p x k Design matrix (p is number of choice alts)} 
  \item{beta}{ k x 1 coefficient vector }
}
\value{
  a list containing:
  \item{y}{ n x 1 vector of multinomial outcomes (1, \ldots, nrow(X)/n)}
  \item{X}{ Design matrix }
  \item{beta}{ coefficient vector }
  \item{prob}{ n x p array of choice probs }
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}
\seealso{ \code{\link{simmnl}} }
\keyword{ models }
