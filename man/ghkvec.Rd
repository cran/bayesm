\name{ghkvec}
\alias{ghkvec}
\title{ Compute GHK approximation to Multivariate Normal Integrals }
\description{
  \code{ghkvec} computes the GHK approximation to the integral of a
  multivariate normal density over a half plane defined by a set
  of truncation points.  
}
\usage{
ghkvec(L, trunpt, above, r)
}
\arguments{
  \item{L}{ lower triangular Cholesky root of Covariance matrix }
  \item{trunpt}{ vector of truncation points}
  \item{above}{ vector of indicators for truncation above(1) or below(0) }
  \item{r}{ number of draws to use in GHK }
}
\value{
  approximation to integral
}
\note{
  \code{ghkvec} can accept a vector of truncations and compute more than one
  integral.  That is, length(trunpt)/length(above) number of different integrals,
  each with the same Sigma and mean 0 but different truncation points. See 
  example below for an example with two integrals at different truncation points.
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 3. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\examples{
##

Sigma=matrix(c(1,.5,.5,1),ncol=2)
L=t(chol(Sigma))
trunpt=c(0,0,1,1)
above=c(1,1)
ghkvec(L,trunpt,above,100)
}
\keyword{ distribution }
