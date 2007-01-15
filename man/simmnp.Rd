\name{simmnp}
\alias{simmnp}
\concept{simulation}
\concept{probit}
\concept{mnp}
\title{ Simulate from Multinomial Probit Model }
\description{
  \code{simmvp} simulates from the multinomial probit model.
}
\usage{
simmnp(X, p, n, beta, sigma)
}
\arguments{
  \item{X}{ n*(p-1) x length(beta) Design matrix }
  \item{p}{ number of choice alternatives }
  \item{n}{ number of observations }
  \item{beta}{ coefficient vector }
  \item{sigma}{ (p-1) x (p-1) covariance matrix }
}
\value{
  a list of 
  \item{y}{n vector of multinomial (1, \ldots,p) outcomes}
  \item{X}{Design matrix}
  \item{beta}{ coefficients}
  \item{sigma}{ covariance matrix }
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 4. \cr
  \url{http://faculty.chicagogsb.edu/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{rmnpGibbs}} }

}
\keyword{ models }
