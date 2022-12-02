\name{simnhlogit}
\alias{simnhlogit}
\concept{logit}
\concept{non-homothetic}

\title{Simulate from Non-homothetic Logit Model}

\description{\code{simnhlogit} simulates from the non-homothetic logit model.}

\usage{simnhlogit(theta, lnprices, Xexpend)}

\arguments{
  \item{theta    }{ coefficient vector }
  \item{lnprices }{ \eqn{n x p} array of prices }
  \item{Xexpend  }{ \eqn{n x k} array of values of expenditure variables}
}

\details{For details on parameterization, see \code{llnhlogit}.}

\value{
A list containing: 
  \item{y        }{\eqn{n x 1} vector of multinomial outcomes (\eqn{1,\ldots,p})}
  \item{Xexpend  }{ expenditure variables}
  \item{lnprices }{ price array }
  \item{theta    }{ coefficients}
  \item{prob     }{ \eqn{n x p} array of choice probabilities}
}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 4, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{llnhlogit}} }

\examples{
N = 1000
p = 3
k = 1

theta = c(rep(1,p), seq(from=-1,to=1,length=p), rep(2,k), 0.5)
lnprices = matrix(runif(N*p), ncol=p)
Xexpend  = matrix(runif(N*k), ncol=k)

simdata = simnhlogit(theta, lnprices, Xexpend)
}

\keyword{models}
