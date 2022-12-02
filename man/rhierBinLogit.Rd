\name{rhierBinLogit}
\alias{rhierBinLogit}
\concept{bayes}
\concept{MCMC}
\concept{hierarchical models}
\concept{binary logit}

\title{MCMC Algorithm for Hierarchical Binary Logit}

\description{
\bold{This function has been deprecated. Please use \code{rhierMnlRwMixture} instead.}

\code{rhierBinLogit} implements an MCMC algorithm for hierarchical binary logits with a normal heterogeneity distribution. This is a hybrid sampler with a RW Metropolis step for unit-level logit parameters. 

\code{rhierBinLogit} is designed for use on choice-based conjoint data with partial profiles. The Design matrix is based on differences of characteristics between two alternatives. See Appendix A of \emph{Bayesian Statistics and Marketing} for details.
}

\usage{rhierBinLogit(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(lgtdata, Z)}
  \item{Prior}{list(Deltabar, ADelta, nu, V)}
  \item{Mcmc }{list(R, keep, sbeta)}
}

\details{
\subsection{Model and Priors}{
  \eqn{y_{hi} = 1} with \eqn{Pr = exp(x_{hi}'\beta_h) / (1+exp(x_{hi}'\beta_h)} and \eqn{\beta_h} is \eqn{nvar x 1} \cr
  \eqn{h = 1, \ldots, length(lgtdata)} units (or "respondents" for survey data)

  \eqn{\beta_h} = ZDelta[h,] + \eqn{u_h} \cr
  Note: here ZDelta refers to \code{Z\%*\%Delta} with ZDelta[h,] the \eqn{h}th row of this product\cr
  Delta is an \eqn{nz x nvar} array

  \eqn{u_h} \eqn{\sim}{~} \eqn{N(0, V_{beta})}.  \cr

  \eqn{delta = vec(Delta)} \eqn{\sim}{~} \eqn{N(vec(Deltabar), V_{beta}(x) ADelta^{-1})}\cr
  \eqn{V_{beta}} \eqn{\sim}{~} \eqn{IW(nu, V)}
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(lgtdata, Z)} [\code{Z} optional]}
  \tabular{ll}{
    \code{lgtdata:       } \tab list of lists with each cross-section unit MNL data \cr
    \code{lgtdata[[h]]$y:} \tab \eqn{n_h x 1} vector of binary outcomes (0,1) \cr
    \code{lgtdata[[h]]$X:} \tab \eqn{n_h x nvar} design matrix for h'th unit \cr
    \code{Z:             } \tab \eqn{nreg x nz} mat of unit chars (def: vector of ones)
  }
  \emph{\code{Prior = list(Deltabar, ADelta, nu, V)} [optional]}
  \tabular{ll}{
    \code{Deltabar:} \tab \eqn{nz x nvar} matrix of prior means (def: 0) \cr
    \code{ADelta:  } \tab prior precision matrix (def: 0.01I) \cr
    \code{nu:      } \tab d.f. parameter for IW prior on normal component Sigma (def: nvar+3) \cr
    \code{V:       } \tab pds location parm for IW prior on normal component Sigma (def: nuI)
  }
  \emph{\code{Mcmc  = list(R, keep, sbeta)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:       } \tab number of MCMC draws \cr
    \code{keep:    } \tab MCMC thinning parm -- keep every \code{keep}th draw (def: 1) \cr
    \code{sbeta:   } \tab scaling parm for RW Metropolis (def: 0.2)
  }
}
}

\value{
  A list containing:
  \item{Deltadraw }{ \eqn{R/keep x nz*nvar} matrix of draws of Delta}
  \item{betadraw  }{ \eqn{nlgt x nvar x R/keep} array of draws of betas}
  \item{Vbetadraw }{ \eqn{R/keep x nvar*nvar} matrix of draws of Vbeta}
  \item{llike     }{ \eqn{R/keep x 1} vector of log-like values}
  \item{reject    }{ \eqn{R/keep x 1} vector of reject rates over nlgt units}
}

\note{Some experimentation with the Metropolis scaling paramter (\code{sbeta}) may be required.}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=10000} else {R=10}
set.seed(66)

nvar = 5              ## number of coefficients
nlgt = 1000           ## number of cross-sectional units
nobs = 10             ## number of observations per unit
nz = 2                ## number of regressors in mixing distribution

Z = matrix(c(rep(1,nlgt),runif(nlgt,min=-1,max=1)), nrow=nlgt, ncol=nz)
Delta = matrix(c(-2, -1, 0, 1, 2, -1, 1, -0.5, 0.5, 0), nrow=nz, ncol=nvar)
iota = matrix(1, nrow=nvar, ncol=1)
Vbeta = diag(nvar) + 0.5*iota\%*\%t(iota)

lgtdata=NULL
for (i in 1:nlgt) { 
  beta = t(Delta)\%*\%Z[i,] + as.vector(t(chol(Vbeta))\%*\%rnorm(nvar))
  X = matrix(runif(nobs*nvar), nrow=nobs, ncol=nvar)
  prob = exp(X\%*\%beta) / (1+exp(X\%*\%beta)) 
  unif = runif(nobs, 0, 1)
  y = ifelse(unif<prob, 1, 0)
  lgtdata[[i]] = list(y=y, X=X, beta=beta)
}

Data1 = list(lgtdata=lgtdata, Z=Z)
Mcmc1 = list(R=R)

out = rhierBinLogit(Data=Data1, Mcmc=Mcmc1)

cat("Summary of Delta draws", fill=TRUE)
summary(out$Deltadraw, tvalues=as.vector(Delta))

cat("Summary of Vbeta draws", fill=TRUE)
summary(out$Vbetadraw, tvalues=as.vector(Vbeta[upper.tri(Vbeta,diag=TRUE)]))

if(0){
## plotting examples
plot(out$Deltadraw,tvalues=as.vector(Delta))
plot(out$betadraw)
plot(out$Vbetadraw,tvalues=as.vector(Vbeta[upper.tri(Vbeta,diag=TRUE)]))
}

}
\seealso{ \code{\link{rhierMnlRwMixture}} }
\keyword{ models}

