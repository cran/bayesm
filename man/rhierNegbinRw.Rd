\name{rhierNegbinRw}
\alias{rhierNegbinRw}
\concept{MCMC}
\concept{hierarchical NBD regression}
\concept{Negative Binomial regression}
\concept{Poisson regression}
\concept{Metropolis algorithm}
\concept{bayes}

\title{MCMC Algorithm for Hierarchical Negative Binomial Regression}

\description{
\code{rhierNegbinRw} implements an MCMC algorithm for the hierarchical Negative Binomial (NBD) regression model. Metropolis steps for each unit-level set of regression parameters are automatically tuned by optimization. Over-dispersion parameter (alpha) is common across units.
}

\usage{rhierNegbinRw(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(regdata, Z)}
  \item{Prior}{list(Deltabar, Adelta, nu, V, a, b)}
  \item{Mcmc }{list(R, keep, nprint, s_beta, s_alpha, alpha, c, Vbeta0, Delta0)}
}

\details{
\subsection{Model and Priors}{
  \eqn{y_i} \eqn{\sim}{~} NBD(mean=\eqn{\lambda}, over-dispersion=alpha) \cr
  \eqn{\lambda = exp(X_i\beta_i)}

  \eqn{\beta_i} \eqn{\sim}{~} \eqn{N(\Delta'z_i,Vbeta)}

  \eqn{vec(\Delta|Vbeta)} \eqn{\sim}{~} \eqn{N(vec(Deltabar), Vbeta (x) Adelta)} \cr
  \eqn{Vbeta}             \eqn{\sim}{~} \eqn{IW(nu, V)} \cr
  \eqn{alpha}             \eqn{\sim}{~} \eqn{Gamma(a, b)} (unless \code{Mcmc$alpha} specified) \cr
  Note: prior mean of \eqn{alpha = a/b}, variance \eqn{= a/(b^2)}
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(regdata, Z)} [\code{Z} optional]}
  \tabular{ll}{
    \code{regdata:        } \tab list of lists with data on each of \code{nreg} units \cr
    \code{regdata[[i]]$X: } \tab \eqn{nobs_i x nvar} matrix of \eqn{X} variables \cr
    \code{regdata[[i]]$y: } \tab \eqn{nobs_i x 1} vector of count responses \cr
    \code{Z:              } \tab \eqn{nreg x nz} matrix of unit characteristics (def: vector of ones)
    }
  \emph{\code{Prior = list(Deltabar, Adelta, nu, V, a, b)} [optional]}
  \tabular{ll}{
    \code{Deltabar:       } \tab \eqn{nz x nvar} prior mean matrix (def: 0) \cr
    \code{Adelta:         } \tab \eqn{nz x nz} PDS prior precision matrix (def: 0.01*I) \cr
    \code{nu:             } \tab d.f. parameter for Inverted Wishart prior (def: nvar+3) \cr
    \code{V:              } \tab location matrix of Inverted Wishart prior (def: nu*I) \cr
    \code{a:              } \tab Gamma prior parameter (def: 0.5) \cr
    \code{b:              } \tab Gamma prior parameter (def: 0.1)
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, s_beta, s_alpha, alpha, c, Vbeta0, Delta0)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:              } \tab number of MCMC draws \cr
    \code{keep:           } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:         } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{s_beta:         } \tab scaling for beta | alpha RW inc cov (def: 2.93/\code{sqrt(nvar)}) \cr
    \code{s_alpha:        } \tab scaling for alpha | beta RW inc cov (def: 2.93) \cr
    \code{alpha:          } \tab over-dispersion parameter (def: alpha ~ Gamma(a,b)) \cr
    \code{c:              } \tab fractional likelihood weighting parm (def: 2) \cr
    \code{Vbeta0:         } \tab starting value for Vbeta (def: I) \cr
    \code{Delta0:         } \tab starting value for Delta (def: 0)
    }
}
}

\value{
  A list containing: 
  \item{llike        }{ \eqn{R/keep x 1} vector of values of log-likelihood}
  \item{betadraw     }{ \eqn{nreg x nvar x R/keep} array of beta draws}
  \item{alphadraw    }{ \eqn{R/keep x 1} vector of alpha draws}
  \item{acceptrbeta  }{ acceptance rate of the beta draws}
  \item{acceptralpha }{ acceptance rate of the alpha draws}
}



\note{
  The NBD regression encompasses Poisson regression in the sense that as alpha goes to infinity the NBD distribution tends to the Poisson.
  
  For "small" values of alpha, the dependent variable can be extremely variable so that a large number of observations may be required to obtain precise inferences. 

  For ease of interpretation, we recommend demeaning \eqn{Z} variables.
}

\seealso{ \code{\link{rnegbinRw}} }

\author{Sridhar Narayanan (Stanford GSB) and Peter Rossi (Anderson School, UCLA), \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

# Simulate from the Negative Binomial Regression
simnegbin = function(X, beta, alpha) {
  lambda = exp(X\%*\%beta)
  y = NULL
  for (j in 1:length(lambda)) {y = c(y, rnbinom(1, mu=lambda[j], size=alpha)) }
  return(y)
  }

nreg = 100        # Number of cross sectional units
T = 50            # Number of observations per unit
nobs = nreg*T
nvar = 2          # Number of X variables
nz = 2            # Number of Z variables
              
## Construct the Z matrix
Z = cbind(rep(1,nreg), rnorm(nreg,mean=1,sd=0.125))

Delta = cbind(c(4,2), c(0.1,-1))
alpha = 5
Vbeta = rbind(c(2,1), c(1,2))

## Construct the regdata (containing X)
simnegbindata = NULL
for (i in 1:nreg) {
    betai = as.vector(Z[i,]\%*\%Delta) + chol(Vbeta)\%*\%rnorm(nvar)
    X = cbind(rep(1,T),rnorm(T,mean=2,sd=0.25))
    simnegbindata[[i]] = list(y=simnegbin(X,betai,alpha), X=X, beta=betai)
}

Beta = NULL
for (i in 1:nreg) {Beta = rbind(Beta,matrix(simnegbindata[[i]]$beta,nrow=1))}
Data1 = list(regdata=simnegbindata, Z=Z)
Mcmc1 = list(R=R)

out = rhierNegbinRw(Data=Data1, Mcmc=Mcmc1)

cat("Summary of Delta draws", fill=TRUE)
summary(out$Deltadraw, tvalues=as.vector(Delta))

cat("Summary of Vbeta draws", fill=TRUE)
summary(out$Vbetadraw, tvalues=as.vector(Vbeta[upper.tri(Vbeta,diag=TRUE)]))

cat("Summary of alpha draws", fill=TRUE)
summary(out$alpha, tvalues=alpha)

## plotting examples
if(0){
  plot(out$betadraw)
  plot(out$alpha,tvalues=alpha)
  plot(out$Deltadraw,tvalues=as.vector(Delta))
}
}

\keyword{models}
