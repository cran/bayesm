\name{rsurGibbs}
\alias{rsurGibbs}
\concept{bayes}
\concept{Gibbs Sampler}
\concept{regression}
\concept{SUR model}
\concept{Seemingly Unrelated Regression}
\concept{MCMC}

\title{Gibbs Sampler for Seemingly Unrelated Regressions (SUR)}

\description{
\code{rsurGibbs} implements a Gibbs Sampler to draw from the posterior of the Seemingly Unrelated Regression (SUR) Model of Zellner.
}

\usage{rsurGibbs(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(regdata)}
  \item{Prior}{list(betabar, A, nu, V)}
  \item{Mcmc }{list(R, keep)}
}

\details{
\subsection{Model and Priors}{
  \eqn{y_i = X_i\beta_i + e_i} with  \eqn{i=1,\ldots,m} for \eqn{m} regressions \cr
  (\eqn{e(k,1), \ldots, e(k,m)})' \eqn{\sim}{~} \eqn{N(0, \Sigma)} with \eqn{k=1, \ldots, n} 

  Can be written as a stacked model: \cr
  \eqn{y = X\beta + e} where \eqn{y} is a \eqn{nobs*m} vector and \eqn{p} = \code{length(beta)} = \code{sum(length(beta_i))}

  Note: must have the same number of observations (\eqn{n}) in each equation but can have a different number of \eqn{X} variables (\eqn{p_i}) for each equation where \eqn{p = \sum p_i}.

  \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar, A^{-1})} \cr
  \eqn{\Sigma} \eqn{\sim}{~} \eqn{IW(nu,V)}
}
\subsection{Argument Details}{
  \emph{\code{Data = list(regdata)}}
    \tabular{ll}{
      \code{regdata: } \tab list of lists, \code{regdata[[i]] = list(y=y_i, X=X_i)}, where \code{y_i} is \eqn{n x 1} and \code{X_i} is \eqn{n x p_i}
    }
  \emph{\code{Prior = list(betabar, A, nu, V)} [optional]}
    \tabular{ll}{
      \code{betabar: } \tab \eqn{p x 1} prior mean (def: 0) \cr
      \code{A:       } \tab \eqn{p x p} prior precision matrix (def: 0.01*I) \cr
      \code{nu:      } \tab d.f. parameter for Inverted Wishart prior (def: m+3) \cr
      \code{V:       } \tab \eqn{m x m} scale parameter for Inverted Wishart prior (def: nu*I)
    }
  \emph{\code{Mcmc = list(R, keep)} [only \code{R} required]}
    \tabular{ll}{
      \code{R:       }\tab number of MCMC draws \cr
      \code{keep:    }\tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
      \code{nprint:  }\tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print)
  }
}
}

\value{
  A list containing: 
  \item{betadraw }{ \eqn{R x p} matrix of betadraws}
  \item{Sigmadraw }{ \eqn{R x (m*m)} array of Sigma draws}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{rmultireg}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=1000} else {R=10}
set.seed(66)

## simulate data from SUR
beta1 = c(1,2)
beta2 = c(1,-1,-2)
nobs = 100
nreg = 2
iota = c(rep(1, nobs))
X1 = cbind(iota, runif(nobs))
X2 = cbind(iota, runif(nobs), runif(nobs))
Sigma = matrix(c(0.5, 0.2, 0.2, 0.5), ncol=2)
U = chol(Sigma)
E = matrix(rnorm(2*nobs),ncol=2)\%*\%U
y1 = X1\%*\%beta1 + E[,1]
y2 = X2\%*\%beta2 + E[,2]

## run Gibbs Sampler
regdata = NULL
regdata[[1]] = list(y=y1, X=X1)
regdata[[2]] = list(y=y2, X=X2)

out = rsurGibbs(Data=list(regdata=regdata), Mcmc=list(R=R))

cat("Summary of beta draws", fill=TRUE)
summary(out$betadraw, tvalues=c(beta1,beta2))

cat("Summary of Sigmadraws", fill=TRUE)
summary(out$Sigmadraw, tvalues=as.vector(Sigma[upper.tri(Sigma,diag=TRUE)]))

## plotting examples
if(0){plot(out$betadraw, tvalues=c(beta1,beta2))}
}

\keyword{regression}
