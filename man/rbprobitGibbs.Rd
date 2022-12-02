\name{rbprobitGibbs}
\alias{rbprobitGibbs}
\concept{bayes}
\concept{MCMC}
\concept{probit}
\concept{Gibbs Sampling}

\title{Gibbs Sampler (Albert and Chib) for Binary Probit}

\description{
\code{rbprobitGibbs} implements the Albert and Chib Gibbs Sampler for the binary probit model.
}

\usage{rbprobitGibbs(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y, X)}
  \item{Prior}{list(betabar, A)}
  \item{Mcmc }{list(R, keep, nprint)}
}

\details{
\subsection{Model and Priors}{
  \eqn{z = X\beta + e} with \eqn{e} \eqn{\sim}{~} \eqn{N(0, I)} \cr
  \eqn{y = 1} if \eqn{z > 0} 

  \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar, A^{-1})}
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(y, X)}}
  \tabular{ll}{
    \code{y:       } \tab \eqn{n x 1} vector of 0/1 outcomes \cr
    \code{X:       } \tab \eqn{n x k} design matrix
    }
  \emph{\code{Prior = list(betabar, A)} [optional]}
  \tabular{ll}{
    \code{betabar: } \tab \eqn{k x 1} prior mean (def: 0) \cr
    \code{A:       } \tab \eqn{k x k} prior precision matrix (def: 0.01*I)
    }
  \emph{\code{Mcmc  = list(R, keep, nprint)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:       } \tab number of MCMC draws \cr
    \code{keep:    } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:  } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print)
    }
}
}

\value{
  A list containing:
  \item{betadraw }{ \eqn{R/keep x k} matrix of betadraws}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{rmnpGibbs}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

## function to simulate from binary probit including x variable
simbprobit = function(X, beta) {
  y = ifelse((X\%*\%beta + rnorm(nrow(X)))<0, 0, 1)
  list(X=X, y=y, beta=beta)
}

nobs = 200
X = cbind(rep(1,nobs), runif(nobs), runif(nobs))
beta = c(0,1,-1)
nvar = ncol(X)
simout = simbprobit(X, beta)

Data1 = list(X=simout$X, y=simout$y)
Mcmc1 = list(R=R, keep=1)

out = rbprobitGibbs(Data=Data1, Mcmc=Mcmc1)
summary(out$betadraw, tvalues=beta)

## plotting example
if(0){plot(out$betadraw, tvalues=beta)}
}

\keyword{models}
