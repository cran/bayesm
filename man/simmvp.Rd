\name{simmvp}
\alias{simmvp}
\concept{probit}
\concept{mvp}
\title{ Simulate from Multivariate Probit Model }
\description{
  \code{simmvp} simulates from the multivariate probit model.
}
\usage{
simmvp(X, p, n, beta, sigma)
}
\arguments{
  \item{X}{ n*p x length(beta) Design matrix }
  \item{p}{ dimension of the MVP}
  \item{n}{ number of observations }
  \item{beta}{ coefficient vector }
  \item{sigma}{ p x p covariance matrix }
}
\value{
  a list of 
  \item{y}{p*n vector of 0/1 binary outcomes}
  \item{X}{Design matrix}
  \item{beta}{ coefficients}
  \item{sigma}{ covariance matrix }
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 4. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{rmvpGibbs}} }

}
\keyword{ models }
