\name{rhierLinearModel}
\alias{rhierLinearModel}
\concept{bayes}
\concept{MCMC}
\concept{Gibbs Sampling}
\concept{hierarchical models}
\concept{linear model}

\title{Gibbs Sampler for Hierarchical Linear Model with Normal Heterogeneity}

\description{
\code{rhierLinearModel} implements a Gibbs Sampler for hierarchical linear models with a normal prior.
}

\usage{rhierLinearModel(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(regdata, Z)}
  \item{Prior}{list(Deltabar, A, nu.e, ssq, nu, V)}
  \item{Mcmc }{list(R, keep, nprint)}
}

\details{
\subsection{Model and Priors}{
  \code{nreg} regression equations with \eqn{nvar} \eqn{X} variables in each equation \cr
  \eqn{y_i = X_i\beta_i + e_i} with \eqn{e_i} \eqn{\sim}{~} \eqn{N(0, \tau_i)} 

  \eqn{\tau_i} \eqn{\sim}{~} nu.e*\eqn{ssq_i/\chi^2_{nu.e}} where \eqn{\tau_i} is the variance of \eqn{e_i}\cr
  \eqn{\beta_i} \eqn{\sim}{~} N(Z\eqn{\Delta}[i,], \eqn{V_{\beta}}) \cr
  Note:  Z\eqn{\Delta} is the matrix \eqn{Z * \Delta} and [i,] refers to \eqn{i}th row of this product\cr

  \eqn{vec(\Delta)} given \eqn{V_{\beta}} \eqn{\sim}{~} \eqn{N(vec(Deltabar), V_{\beta}(x) A^{-1})}\cr
  \eqn{V_{\beta}} \eqn{\sim}{~} \eqn{IW(nu,V)} \cr
  \eqn{Delta, Deltabar} are \eqn{nz x nvar}; \eqn{A} is \eqn{nz x nz}; \eqn{V_{\beta}} is \eqn{nvar x nvar}.
  
  Note: if you don't have any \eqn{Z} variables, omit \eqn{Z} in the \code{Data} argument and 
  a vector of ones will be inserted; the matrix \eqn{\Delta} will be \eqn{1 x nvar} and should 
  be interpreted as the mean of all unit \eqn{\beta}s.
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(regdata, Z)} [\code{Z} optional]}
  \tabular{ll}{
    \code{regdata:        } \tab list of lists with \eqn{X} and \eqn{y} matrices for each of \code{nreg=length(regdata)} regressions \cr
    \code{regdata[[i]]$X: } \tab \eqn{n_i x nvar} design matrix for equation \eqn{i} \cr
    \code{regdata[[i]]$y: } \tab \eqn{n_i x 1} vector of observations for equation \eqn{i} \cr
    \code{Z:              } \tab \eqn{nreg x nz} matrix of unit characteristics (def: vector of ones)
    }
  \emph{\code{Prior = list(Deltabar, A, nu.e, ssq, nu, V)}  [optional]}
  \tabular{ll}{
    \code{Deltabar:       } \tab \eqn{nz x nvar} matrix of prior means (def: 0) \cr
    \code{A:              } \tab \eqn{nz x nz} matrix for prior precision (def: 0.01I) \cr
    \code{nu.e:           } \tab d.f. parameter for regression error variance prior (def: 3) \cr
    \code{ssq:            } \tab scale parameter for regression error var prior (def: \code{var(y_i)}) \cr
    \code{nu:             } \tab d.f. parameter for Vbeta prior (def: nvar+3) \cr
    \code{V:              } \tab Scale location matrix for Vbeta prior (def: nu*I)
    }
  \emph{\code{Mcmc  = list(R, keep, nprint)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:              } \tab number of MCMC draws \cr
    \code{keep:           } \tab MCMC thinning parm -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:         } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print)
    }
}
}

\value{
  A list containing:
  \item{betadraw  }{\eqn{nreg x nvar x R/keep} array of individual regression coef draws}
  \item{taudraw   }{\eqn{R/keep x nreg} matrix of error variance draws}
  \item{Deltadraw }{\eqn{R/keep x nz*nvar} matrix of Deltadraws}
  \item{Vbetadraw }{\eqn{R/keep x nvar*nvar} matrix of Vbeta draws}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{rhierLinearMixture}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

nreg = 100
nobs = 100
nvar = 3
Vbeta = matrix(c(1, 0.5, 0, 0.5, 2, 0.7, 0, 0.7, 1), ncol=3)

Z = cbind(c(rep(1,nreg)), 3*runif(nreg))
Z[,2] = Z[,2] - mean(Z[,2])
nz = ncol(Z)
Delta = matrix(c(1,-1,2,0,1,0), ncol=2)
Delta = t(Delta) # first row of Delta is means of betas
Beta = matrix(rnorm(nreg*nvar),nrow=nreg)\%*\%chol(Vbeta) + Z\%*\%Delta

tau = 0.1
iota = c(rep(1,nobs))
regdata = NULL
for (reg in 1:nreg) { 
  X = cbind(iota, matrix(runif(nobs*(nvar-1)),ncol=(nvar-1)))
	y = X\%*\%Beta[reg,] + sqrt(tau)*rnorm(nobs)
	regdata[[reg]] = list(y=y, X=X) 
}

Data1 = list(regdata=regdata, Z=Z)
Mcmc1 = list(R=R, keep=1)

out = rhierLinearModel(Data=Data1, Mcmc=Mcmc1)

cat("Summary of Delta draws", fill=TRUE)
summary(out$Deltadraw, tvalues=as.vector(Delta))

cat("Summary of Vbeta draws", fill=TRUE)
summary(out$Vbetadraw, tvalues=as.vector(Vbeta[upper.tri(Vbeta,diag=TRUE)]))

## plotting examples
if(0){
  plot(out$betadraw)
  plot(out$Deltadraw)
}
}

\keyword{regression}
