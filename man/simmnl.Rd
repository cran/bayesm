\name{simmnl}
\alias{simmnl}
\concept{logit}
\concept{mnl}
\title{ Simulate from Multinomial Logit Model }
\description{
  \code{simmnl} simulates from the MNL model.
}
\usage{
simmnl(p, n, beta)
}
\arguments{
  \item{p}{ number choice alternatives }
  \item{n}{ number of observations }
  \item{beta}{ MNL coefficient vector }
}
\details{
  \code{simmnl} will simulate two uniformly distributed X vars and add intercepts.
}

\value{
  \item{y}{ n x 1 vector of multinomial outcomes (1, \ldots, p)}
  \item{X} { n x k design matrix }
  \item{beta}{beta vector }
  \item{prob}{n x j array of choice probabilities }
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby, and McCulloch. \cr
  \url{http://faculty.chicagogsb.edu/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}
\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{llmnl}}, \code{\link{rmnlIndepMetrop}} }
}
\keyword{ models }
