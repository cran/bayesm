\name{runireg}
\alias{runireg}
\concept{bayes}
\concept{regression}

\title{IID Sampler for Univariate Regression}

\description{
\code{runireg} implements an iid sampler to draw from the posterior of a univariate regression with a conjugate prior.
}

\usage{runireg(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y, X)}
  \item{Prior}{list(betabar, A, nu, ssq)}
  \item{Mcmc }{list(R, keep, nprint)}
}

\details{
\subsection{Model and Priors}{
  \eqn{y = X\beta + e} with \eqn{e} \eqn{\sim}{~} \eqn{N(0, \sigma^2)}

  \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar, \sigma^2*A^{-1})} \cr
  \eqn{\sigma^2} \eqn{\sim}{~} \eqn{(nu*ssq)/\chi^2_{nu}}
}
\subsection{Argument Details}{
  \emph{\code{Data = list(y, X)}}
    \tabular{ll}{
      \code{y:       } \tab \eqn{n x 1} vector of observations \cr
      \code{X:       } \tab \eqn{n x k} design matrix
    }
  \emph{\code{Prior = list(betabar, A, nu, ssq)} [optional]}
    \tabular{ll}{
      \code{betabar: } \tab \eqn{k x 1} prior mean (def: 0) \cr
      \code{A:       } \tab \eqn{k x k} prior precision matrix (def: 0.01*I) \cr
      \code{nu:      } \tab d.f. parameter for Inverted Chi-square prior (def: 3) \cr
      \code{ssq:     } \tab scale parameter for Inverted Chi-square prior (def: \code{var(y)})
    }
  \emph{\code{Mcmc = list(R, keep, nprint)} [only \code{R} required]}
    \tabular{ll}{
      \code{R:       } \tab number of draws \cr
      \code{keep:    } \tab thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
      \code{nprint:  } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print)
    }
}
}

\value{
  A list containing: 
  \item{betadraw }{ \eqn{R x k} matrix of betadraws }
  \item{sigmasqdraw }{ \eqn{R x 1} vector of sigma-sq draws}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 2, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{runiregGibbs}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

n = 200
X = cbind(rep(1,n), runif(n))
beta = c(1,2)
sigsq = 0.25
y = X\%*\%beta + rnorm(n,sd=sqrt(sigsq))

out = runireg(Data=list(y=y,X=X), Mcmc=list(R=R))

cat("Summary of beta and Sigmasq draws", fill=TRUE)
summary(out$betadraw, tvalues=beta)
summary(out$sigmasqdraw, tvalues=sigsq)

## plotting examples
if(0){plot(out$betadraw)}
}

\keyword{regression}
