\name{init.rmultiregfp}
\alias{init.rmultiregfp}

\title{ Initialize Variables for Multivariate Regression Draw }
\description{
  \code{init.rmultiregfp} initializes variables which can be pre-computed
  for draws from the posterior of a multivariate regression model.  
  \code{init.rmultiregfp} should be called prior to use of \code{rmultiregfp}
}
\usage{
init.rmultiregfp(X, A, Bbar, nu, V)
}
\arguments{
  \item{X}{ Design matrix }
  \item{A}{ Prior Precision matrix (m x m) }
  \item{Bbar}{ Prior mean matrix (m x k) }
  \item{nu}{ degrees of freedom parmeter for Sigma prior }
  \item{V}{ location parameter for Sigma prior }
}
\details{
  model: \eqn{Y = XB + U}. \eqn{u_i} \eqn{\sim}{~} \eqn{N(0,\Sigma)}.  \eqn{u_i} is the ith row of U.  \eqn{Y} is n x m. \eqn{X} is n x k. \eqn{B} is k x m. 

  priors: \eqn{vec(B)} \eqn{\sim}{~} \eqn{N(vec(Bbar,\Sigma (x) A^{-1})} \cr
          \eqn{\Sigma} \eqn{\sim}{~} \eqn{IW(nu,V)}.
}
\value{
  A list containing \ldots
  \item{IR }{Inverse of Cholesky Root of \eqn{(X'X + A)}}
  \item{RA }{ Cholesky root of A }
  \item{RABbar}{RA \%*\% Bbar}
  \item{nu}{d.f. parm for IWishart prior}
  \item{V}{location matrix for IWishart prior}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 2. \cr
  \url{http://faculty.chicagogsb.edu/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\seealso{ \code{\link{rmultiregfp}}}

\keyword{ utilities }
\keyword{ regression }
