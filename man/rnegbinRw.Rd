\name{rnegbinRw}
\alias{rnegbinRw}
\concept{MCMC}
\concept{NBD regression}
\concept{Negative Binomial regression}
\concept{Poisson regression}
\concept{Metropolis algorithm}
\concept{bayes}

\title{MCMC Algorithm for Negative Binomial Regression}

\description{
\code{rnegbinRw} implements a Random Walk Metropolis Algorithm for the Negative Binomial (NBD) regression model where \eqn{\beta|\alpha} and \eqn{\alpha|\beta} are drawn with two different random walks.
}

\usage{rnegbinRw(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y, X)}
  \item{Prior}{list(betabar, A, a, b)}
  \item{Mcmc }{list(R, keep, nprint, s_beta, s_alpha, beta0, alpha)}
}

\details{
\subsection{Model and Priors}{
  \eqn{y} \eqn{\sim}{~} \eqn{NBD(mean=\lambda, over-dispersion=alpha)}  \cr
  \eqn{\lambda = exp(x'\beta)}
  
  \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar, A^{-1})} \cr
  \eqn{alpha} \eqn{\sim}{~} \eqn{Gamma(a, b)} (unless \code{Mcmc$alpha} specified) \cr
  Note: prior mean of \eqn{alpha = a/b}, \eqn{variance = a/(b^2)}
}
\subsection{Argument Details}{
  \emph{\code{Data = list(y, X)}}
  \tabular{ll}{
    \code{y:       } \tab \eqn{n x 1} vector of counts (\eqn{0,1,2,\ldots}) \cr
    \code{X:       } \tab \eqn{n x k} design matrix
    }
  \emph{\code{Prior = list(betabar, A, a, b)} [optional]}
  \tabular{ll}{
    \code{betabar: } \tab \eqn{k x 1} prior mean (def: 0) \cr
    \code{A:       } \tab \eqn{k x k} PDS prior precision matrix (def: 0.01*I) \cr
    \code{a:       } \tab Gamma prior parameter (not used if \code{Mcmc$alpha} specified) (def: 0.5) \cr
    \code{b:       } \tab Gamma prior parameter (not used if \code{Mcmc$alpha} specified) (def: 0.1)
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, s_beta, s_alpha, beta0, alpha)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:       } \tab number of MCMC draws \cr
    \code{keep:    } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:  } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{s_beta:  } \tab scaling for beta | alpha RW inc cov matrix (def: 2.93/\code{sqrt(k)}) \cr
    \code{s_alpha: } \tab scaling for alpha | beta RW inc cov matrix (def: 2.93) \cr
    \code{alpha:   } \tab over-dispersion parameter (def: alpha ~ Gamma(a,b))
    }
}
}

\value{
  A list containing: 
  \item{betadraw     }{\eqn{R/keep x k} matrix of beta draws}
  \item{alphadraw    }{\eqn{R/keep x 1} vector of alpha draws}
  \item{llike        }{\eqn{R/keep x 1} vector of log-likelihood values evaluated at each draw}
  \item{acceptrbeta  }{acceptance rate of the beta draws}
  \item{acceptralpha }{acceptance rate of the alpha draws}
}

\note{
The NBD regression encompasses Poisson regression in the sense that as alpha goes to infinity the NBD distribution tends toward the Poisson. For "small" values of alpha, the dependent variable can be extremely variable so that a large number of observations may be required to obtain precise inferences.
}

\author{Sridhar Narayanan (Stanford GSB) and Peter Rossi (Anderson School, UCLA), \email{perossichi@gmail.com}.}

\references{For further discussion, see \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, McCulloch.}

\seealso{ \code{\link{rhierNegbinRw}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0)  {R=1000} else {R=10}
set.seed(66)

simnegbin = function(X, beta, alpha) {
  # Simulate from the Negative Binomial Regression
  lambda = exp(X\%*\%beta)
  y = NULL
  for (j in 1:length(lambda)) { y = c(y, rnbinom(1, mu=lambda[j], size=alpha)) }
  return(y)
}

nobs = 500
nvar = 2 # Number of X variables
alpha = 5
Vbeta = diag(nvar)*0.01

# Construct the regdata (containing X)
simnegbindata = NULL
beta = c(0.6, 0.2)
X = cbind(rep(1,nobs), rnorm(nobs,mean=2,sd=0.5))
simnegbindata = list(y=simnegbin(X,beta,alpha), X=X, beta=beta)

Data1 = simnegbindata
Mcmc1 = list(R=R)

out = rnegbinRw(Data=Data1, Mcmc=list(R=R))

cat("Summary of alpha/beta draw", fill=TRUE)
summary(out$alphadraw, tvalues=alpha)
summary(out$betadraw, tvalues=beta)

## plotting examples
if(0){plot(out$betadraw)}
}

\keyword{models}

