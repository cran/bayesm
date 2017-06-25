\name{rhierLinearMixture}
\alias{rhierLinearMixture}
\concept{bayes}
\concept{MCMC}
\concept{Gibbs Sampling}
\concept{mixture of normals}
\concept{normal mixture}
\concept{heterogeneity}
\concept{regresssion}
\concept{hierarchical models}
\concept{linear model}

\title{Gibbs Sampler for Hierarchical Linear Model with Mixture-of-Normals Heterogeneity}

\description{
\code{rhierLinearMixture} implements a Gibbs Sampler for hierarchical linear models with a mixture-of-normals prior.
}

\usage{rhierLinearMixture(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(regdata, Z)}
  \item{Prior}{list(deltabar, Ad, mubar, Amu, nu, V, nu.e, ssq, ncomp)}
  \item{Mcmc }{list(R, keep, nprint)}
}

\details{
\subsection{Model and Priors}{
  \code{nreg} regression equations with \code{nvar} as the number of \eqn{X} vars in each equation \cr
  \eqn{y_i = X_i\beta_i + e_i} with \eqn{e_i} \eqn{\sim}{~} \eqn{N(0, \tau_i)}  

  \eqn{\tau_i} \eqn{\sim}{~} \eqn{nu.e*ssq_i/\chi^2_{nu.e}} where  \eqn{\tau_i} is the variance of \eqn{e_i}\cr
  \eqn{B = Z\Delta + U} or \eqn{\beta_i = \Delta' Z[i,]' + u_i} \cr
  \eqn{\Delta} is an \eqn{nz x nvar} matrix \cr

  \eqn{Z} should \emph{not} include an intercept and should be centered for ease of interpretation. 
  The mean of each of the \code{nreg} \eqn{\beta}s is the mean of the normal mixture. 
  Use \code{summary()} to compute this mean from the \code{compdraw} output.

  \eqn{u_i} \eqn{\sim}{~} \eqn{N(\mu_{ind}, \Sigma_{ind})}\cr
  \eqn{ind} \eqn{\sim}{~} \eqn{multinomial(pvec)} \cr

  \eqn{pvec} \eqn{\sim}{~} \eqn{dirichlet(a)}\cr
  \eqn{delta = vec(\Delta)} \eqn{\sim}{~} \eqn{N(deltabar, A_d^{-1})}\cr
  \eqn{\mu_j} \eqn{\sim}{~} \eqn{N(mubar, \Sigma_j(x) Amu^{-1})}\cr
  \eqn{\Sigma_j} \eqn{\sim}{~} \eqn{IW(nu, V)} \cr
  
  Be careful in assessing the prior parameter \code{Amu}: 0.01 can be too small for some applications. 
  See chapter 5 of Rossi et al for full discussion.\cr
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(regdata, Z)} [\code{Z} optional]}
  \tabular{ll}{
    \code{regdata:        } \tab list of lists with \eqn{X} and \eqn{y} matrices for each of \code{nreg=length(regdata)} regressions \cr
    \code{regdata[[i]]$X: } \tab \eqn{n_i x nvar} design matrix for equation \eqn{i}  \cr
    \code{regdata[[i]]$y: } \tab \eqn{n_i x 1} vector of observations for equation \eqn{i} \cr
    \code{Z:              } \tab \eqn{nreg x nz} matrix of unit characteristics (def: vector of ones)
    }
  \emph{\code{Prior = list(deltabar, Ad, mubar, Amu, nu, V, nu.e, ssq, ncomp)} [all but \code{ncomp} are optional]}
  \tabular{ll}{
    \code{deltabar:       } \tab \eqn{nz x nvar} vector of prior means (def: 0) \cr
    \code{Ad:             } \tab prior precision matrix for vec(Delta) (def: 0.01*I) \cr
    \code{mubar:          } \tab \eqn{nvar x 1} prior mean vector for normal component mean (def: 0) \cr
    \code{Amu:            } \tab prior precision for normal component mean (def: 0.01) \cr
    \code{nu:             } \tab d.f. parameter for IW prior on normal component Sigma (def: nvar+3) \cr
    \code{V:              } \tab PDS location parameter for IW prior on normal component Sigma (def: nu*I) \cr
    \code{nu.e:           } \tab d.f. parameter for regression error variance prior (def: 3) \cr
    \code{ssq:            } \tab scale parameter for regression error variance prior (def: \code{var(y_i)}) \cr
    \code{a:              } \tab Dirichlet prior parameter (def: 5) \cr
    \code{ncomp:          } \tab number of components used in normal mixture
    }
  \emph{\code{Mcmc  = list(R, keep, nprint)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:              } \tab number of MCMC draws \cr
    \code{keep:           } \tab MCMC thinning parm -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:         } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print)
    }
}
\subsection{\code{nmix} Details}{
  \code{nmix} is a list with 3 components. Several functions in the \code{bayesm} package that involve a Dirichlet Process or mixture-of-normals return \code{nmix}. Across these functions, a common structure is used for \code{nmix} in order to utilize generic summary and plotting functions. 
  \tabular{ll}{
  \code{probdraw:} \tab \eqn{ncomp x R/keep} matrix that reports the probability that each draw came from a particular component \cr
  \code{zdraw:   } \tab \eqn{R/keep x nobs} matrix that indicates which component each draw is assigned to (here, null) \cr
  \code{compdraw:} \tab A list of \eqn{R/keep} lists of \eqn{ncomp} lists. Each of the inner-most lists has 2 elemens: a vector of draws for \code{mu} and a matrix of draws for the Cholesky root of \code{Sigma}.
  }
}
}

\value{
  A list containing:
  \item{taudraw   }{\eqn{R/keep x nreg} matrix of error variance draws}
  \item{betadraw  }{\eqn{nreg x nvar x R/keep} array of individual regression coef draws}
  \item{Deltadraw }{\eqn{R/keep x nz*nvar} matrix of Deltadraws}
  \item{nmix      }{a list containing: \code{probdraw}, \code{zdraw}, \code{compdraw} (see \dQuote{\code{nmix} Details} section)}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. \cr \url{http://www.perossi.org/home/bsm-1}}

\seealso{ \code{\link{rhierLinearModel}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

nreg = 300
nobs = 500
nvar = 3
nz = 2

Z = matrix(runif(nreg*nz), ncol=nz) 
Z = t(t(Z) - apply(Z,2,mean))
Delta = matrix(c(1,-1,2,0,1,0), ncol=nz)
tau0 = 0.1
iota = c(rep(1,nobs))

## create arguments for rmixture

tcomps = NULL
a = matrix(c(1,0,0,0.5773503,1.1547005,0,-0.4082483,0.4082483,1.2247449), ncol=3)
tcomps[[1]] = list(mu=c(0,-1,-2),   rooti=a) 
tcomps[[2]] = list(mu=c(0,-1,-2)*2, rooti=a)
tcomps[[3]] = list(mu=c(0,-1,-2)*4, rooti=a)
tpvec = c(0.4, 0.2, 0.4)

## simulated data with Z
regdata = NULL
betas = matrix(double(nreg*nvar), ncol=nvar)
tind = double(nreg)

for (reg in 1:nreg) {
  tempout = rmixture(1,tpvec,tcomps)
  betas[reg,] = Delta\%*\%Z[reg,] + as.vector(tempout$x)
  tind[reg] = tempout$z
  X = cbind(iota, matrix(runif(nobs*(nvar-1)),ncol=(nvar-1)))
  tau = tau0*runif(1,min=0.5,max=1)
  y = X\%*\%betas[reg,] + sqrt(tau)*rnorm(nobs)
  regdata[[reg]] = list(y=y, X=X, beta=betas[reg,], tau=tau)
}

## run rhierLinearMixture

Data1 = list(regdata=regdata, Z=Z)
Prior1 = list(ncomp=3)
Mcmc1 = list(R=R, keep=1)

out1 = rhierLinearMixture(Data=Data1, Prior=Prior1, Mcmc=Mcmc1)

cat("Summary of Delta draws", fill=TRUE)
summary(out1$Deltadraw, tvalues=as.vector(Delta))

cat("Summary of Normal Mixture Distribution", fill=TRUE)
summary(out1$nmix)

## plotting examples
if(0){
  plot(out1$betadraw)
  plot(out1$nmix)
  plot(out1$Deltadraw)
}
}

\keyword{ regression }
