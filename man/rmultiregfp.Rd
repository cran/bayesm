\name{rmultiregfp}
\alias{rmultiregfp}

\title{ Draw from the Posterior of a Multivariate Regression }
\description{
  \code{ rmultiregfp} draws from the posterior of a Multivariate Regression model with a natural conjugate
  prior.
}
\usage{
rmultiregfp(Y, X, Fparm)
}

\arguments{
  \item{Y}{ n x m matrix of observations on m dep vars }
  \item{X}{ n x k matrix of observations on indep vars (supply intercept) }
  \item{Fparm}{ a list of prior parameters prepared by \code{init.rmultiregfp}}
}
\details{
  Model: \eqn{Y=XB+U}.  \eqn{cov(u_i) = Sigma}.  \eqn{B} is k x m matrix of coefficients. 
         \eqn{Sigma} is an m x m covariance matrix.

  Priors:  \eqn{beta} given \eqn{Sigma} \eqn{\sim}{~} \eqn{N(betabar,Sigma (x) A^{-1})}. 
  \eqn{betabar=vec(Bbar)};  \eqn{beta = vec(B)}. \cr
  \eqn{Sigma} \eqn{\sim}{~} IW(nu,V). 

  prepare Fparm by call \code{\link{init.rmultiregfp}}
}
\value{
  A list of the components of a draw from the posterior
  \item{B }{ draw of regression coefficient matrix }
  \item{Sigma }{ draw of \eqn{Sigma} }
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

\seealso{ \code{\link{rmultireg}},\code{\link{init.rmultiregfp}}}
\keyword{ regression }
