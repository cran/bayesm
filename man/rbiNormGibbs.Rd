\name{rbiNormGibbs}
\alias{rbiNormGibbs}
\concept{bayes}
\concept{Gibbs Sampling}
\concept{MCMC}
\concept{normal distribution}

\title{Illustrate Bivariate Normal Gibbs Sampler}

\description{
\code{rbiNormGibbs} implements a Gibbs Sampler for the bivariate normal distribution. Intermediate moves are plotted and the output is contrasted with the iid sampler. This function is designed for illustrative/teaching purposes.}

\usage{rbiNormGibbs(initx=2, inity=-2, rho, burnin=100, R=500)}

\arguments{
  \item{initx  }{ initial value of parameter on x axis (def: 2)}
  \item{inity  }{ initial value of parameter on y axis (def: -2)}
  \item{rho    }{ correlation for bivariate normals}
  \item{burnin }{ burn-in number of draws (def: 100)}
  \item{R      }{ number of MCMC draws (def: 500)}
}

\details{
  \eqn{(\theta_1, \theta_2) ~ N( (0,0), \Sigma}) with \eqn{\Sigma} = \code{matrix(c(1,rho,rho,1),ncol=2)}
}

\value{
 \eqn{R x 2} matrix of draws
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapters 2 and 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
\dontrun{out=rbiNormGibbs(rho=0.95)}
}

\keyword{distribution}
