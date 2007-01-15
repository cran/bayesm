\name{simnhlogit}
\alias{simnhlogit}
\concept{logit}
\concept{non-homothetic}
\title{ Simulate from Non-homothetic Logit Model }
\description{
  \code{simnhlogit} simulates from the non-homothetic logit model
}
\usage{
simnhlogit(theta, lnprices, Xexpend)
}
\arguments{
  \item{theta}{ coefficient vector }
  \item{lnprices}{ n x p array of prices }
  \item{Xexpend}{ n x k array of values of expenditure variables}
}
\details{
  For detail on parameterization, see \code{llnhlogit}.
}
\value{
  a list containing: 
  \item{y}{n x 1 vector of multinomial outcomes (1, \ldots, p)}
  \item{Xexpend}{expenditure variables}
  \item{lnprices}{ price array }
  \item{theta}{coefficients}
  \item{prob}{n x p array of choice probabilities}
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

\seealso{ \code{\link{llnhlogit}} }

\keyword{ models }
