\name{rhierMnlRwMixture}
\alias{rhierMnlRwMixture}
\concept{bayes}
\concept{MCMC}
\concept{Multinomial Logit}
\concept{mixture of normals}
\concept{normal mixture}
\concept{heterogeneity}
\concept{hierarchical models}

\title{MCMC Algorithm for Hierarchical Multinomial Logit with Mixture-of-Normals Heterogeneity}

\description{
\code{rhierMnlRwMixture} is a MCMC algorithm for a hierarchical multinomial logit with a mixture of normals heterogeneity distribution. This is a hybrid Gibbs Sampler with a RW Metropolis step for the MNL coefficients for each panel unit. 
}

\usage{rhierMnlRwMixture(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(lgtdata, Z, p)}
  \item{Prior}{list(a, deltabar, Ad, mubar, Amu, nu, V, a, ncomp, SignRes)}
  \item{Mcmc }{list(R, keep, nprint, s, w)}
}

\details{
\subsection{Model and Priors}{
  \eqn{y_i} \eqn{\sim}{~} \eqn{MNL(X_i,\beta_i)} with \eqn{i = 1, \ldots,} length(lgtdata) 
  and where \eqn{\beta_i} is \eqn{nvar x 1}

  \eqn{\beta_i} = \eqn{Z\Delta}[i,] + \eqn{u_i} \cr
  Note:  Z\eqn{\Delta} is the matrix Z * \eqn{\Delta} and [i,] refers to \eqn{i}th row of this product \cr
  Delta is an \eqn{nz x nvar} array

  \eqn{u_i} \eqn{\sim}{~} \eqn{N(\mu_{ind},\Sigma_{ind})} with \eqn{ind} \eqn{\sim}{~} multinomial(pvec)

  \eqn{pvec}                \eqn{\sim}{~} dirichlet(a) \cr
  \eqn{delta = vec(\Delta)} \eqn{\sim}{~} \eqn{N(deltabar, A_d^{-1})} \cr
  \eqn{\mu_j}               \eqn{\sim}{~} \eqn{N(mubar, \Sigma_j (x) Amu^{-1})} \cr
  \eqn{\Sigma_j}            \eqn{\sim}{~} \eqn{IW(nu, V)}
  
  Note: \eqn{Z} should NOT include an intercept and is centered for ease of interpretation. 
  The mean of each of the \code{nlgt} \eqn{\beta}s is the mean of the normal mixture. 
  Use \code{summary()} to compute this mean from the \code{compdraw} output.\cr
  
  Be careful in assessing prior parameter \code{Amu}: 0.01 is too small for many applications. 
  See chapter 5 of Rossi et al for full discussion.
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(lgtdata, Z, p)} [\code{Z} optional]}
  \tabular{ll}{
    \code{lgtdata:        } \tab list of \code{nlgt=length(lgtdata)} lists with each cross-section unit MNL data \cr
    \code{lgtdata[[i]]$y: } \tab \eqn{n_i x 1} vector of multinomial outcomes (1, \ldots, p) \cr
    \code{lgtdata[[i]]$X: } \tab \eqn{n_i*p x nvar} design matrix for \eqn{i}th unit \cr
    \code{Z:              } \tab \eqn{nreg x nz} matrix of unit chars (def: vector of ones) \cr
    \code{p:              } \tab number of choice alternatives 
    }
  \emph{\code{Prior = list(a, deltabar, Ad, mubar, Amu, nu, V, a, ncomp, SignRes)} [all but \code{ncomp} are optional]}
  \tabular{ll}{
    \code{a:              } \tab \eqn{ncomp x 1} vector of Dirichlet prior parameters (def: \code{rep(5,ncomp)}) \cr
    \code{deltabar:       } \tab \eqn{nz*nvar x 1} vector of prior means (def: 0) \cr
    \code{Ad:             } \tab prior precision matrix for vec(D) (def: 0.01*I) \cr
    \code{mubar:          } \tab \eqn{nvar x 1} prior mean vector for normal component mean (def: 0 if unrestricted; 2 if restricted) \cr
    \code{Amu:            } \tab prior precision for normal component mean (def: 0.01 if unrestricted; 0.1 if restricted) \cr
    \code{nu:             } \tab d.f. parameter for IW prior on normal component Sigma (def: nvar+3 if unrestricted; nvar+15 if restricted) \cr
    \code{V:              } \tab PDS location parameter for IW prior on normal component Sigma (def: nu*I if unrestricted; nu*D if restricted with d_pp = 4 if unrestricted and d_pp = 0.01 if restricted) \cr
    \code{ncomp:          } \tab number of components used in normal mixture \cr
    \code{SignRes:        } \tab \eqn{nvar x 1} vector of sign restrictions on the coefficient estimates (def: \code{rep(0,nvar)})
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, s, w)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:              } \tab number of MCMC draws \cr
    \code{keep:           } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:         } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{s:              } \tab scaling parameter for RW Metropolis (def: 2.93/\code{sqrt(nvar)}) \cr
    \code{w:              } \tab fractional likelihood weighting parameter (def: 0.1)
    }
}
\subsection{Sign Restrictions}{
  If \eqn{\beta_ik} has a sign restriction: \eqn{\beta_ik = SignRes[k] * exp(\beta*_ik)}
  
  To use sign restrictions on the coefficients, \code{SignRes} must be an \eqn{nvar x 1} vector containing values of either 0, -1, or 1. The value 0 means there is no sign restriction, -1 ensures that the coefficient is \emph{negative}, and 1 ensures that the coefficient is \emph{positive}. For example, if \code{SignRes = c(0,1,-1)}, the first coefficient is unconstrained, the second will be positive, and the third will be negative.
  
  The sign restriction is implemented such that if the the \eqn{k}'th \eqn{\beta} has a non-zero sign restriction (i.e., it is constrained), we have \eqn{\beta_k = SignRes[k] * exp(\beta*_k)}.
  
  The sign restrictions (if used) will be reflected in the \code{betadraw} output. However, the unconstrained mixture components are available in \code{nmix}. \bold{Important:} Note that draws from \code{nmix} are distributed according to the mixture of normals but \bold{not} the coefficients in \code{betadraw}.
  
  Care should be taken when selecting priors on any sign restricted coefficients. See related vignette for additional information.
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
  \item{Deltadraw }{ \eqn{R/keep x nz*nvar} matrix of draws of Delta, first row is initial value}
  \item{betadraw  }{ \eqn{nlgt x nvar x R/keep} array of beta draws}
  \item{nmix      }{ a list containing: \code{probdraw}, \code{zdraw}, \code{compdraw} (see \dQuote{\code{nmix} Details} section)}
  \item{loglike   }{ \eqn{R/keep x 1} vector of log-likelihood for each kept draw}
  \item{SignRes   }{ \eqn{nvar x 1} vector of sign restrictions}
}

\note{
  Note: as of version 2.0-2 of \code{bayesm}, the fractional weight parameter has been changed to a weight between 0 and 1. 
  \eqn{w} is the fractional weight on the normalized pooled likelihood. This differs from what is in Rossi et al chapter 5, i.e.

  \eqn{like_i^{(1-w)} x like_pooled^{((n_i/N)*w)}}

  Large \code{R} values may be required (>20,000).
} 

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. \cr \url{http://www.perossi.org/home/bsm-1}}
  
\seealso{ \code{\link{rmnlIndepMetrop}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=10000} else {R=10}
set.seed(66)

p = 3                                # num of choice alterns
ncoef = 3  
nlgt = 300                           # num of cross sectional units
nz = 2
Z = matrix(runif(nz*nlgt),ncol=nz)
Z = t(t(Z) - apply(Z,2,mean))        # demean Z
ncomp = 3                            # num of mixture components
Delta = matrix(c(1,0,1,0,1,2),ncol=2)

comps=NULL
comps[[1]] = list(mu=c(0,-1,-2),   rooti=diag(rep(1,3)))
comps[[2]] = list(mu=c(0,-1,-2)*2, rooti=diag(rep(1,3)))
comps[[3]] = list(mu=c(0,-1,-2)*4, rooti=diag(rep(1,3)))
pvec = c(0.4, 0.2, 0.4)

##  simulate from MNL model conditional on X matrix
simmnlwX= function(n,X,beta) {
  k = length(beta)
  Xbeta = X\%*\%beta
  j = nrow(Xbeta) / n
  Xbeta = matrix(Xbeta, byrow=TRUE, ncol=j)
  Prob = exp(Xbeta)
  iota = c(rep(1,j))
  denom = Prob\%*\%iota
  Prob = Prob/as.vector(denom)
  y = vector("double",n)
  ind = 1:j
  for (i in 1:n) { 
    yvec = rmultinom(1, 1, Prob[i,])
    y[i] = ind\%*\%yvec
  }
  return(list(y=y, X=X, beta=beta, prob=Prob))
}

## simulate data
simlgtdata = NULL
ni = rep(50, 300)
for (i in 1:nlgt) {
  betai = Delta\%*\%Z[i,] + as.vector(rmixture(1,pvec,comps)$x)
   Xa = matrix(runif(ni[i]*p,min=-1.5,max=0), ncol=p)
   X = createX(p, na=1, nd=NULL, Xa=Xa, Xd=NULL, base=1)
   outa = simmnlwX(ni[i], X, betai)
   simlgtdata[[i]] = list(y=outa$y, X=X, beta=betai)
}

## plot betas
if(0){
  bmat = matrix(0, nlgt, ncoef)
  for(i in 1:nlgt) {bmat[i,] = simlgtdata[[i]]$beta}
  par(mfrow = c(ncoef,1))
  for(i in 1:ncoef) { hist(bmat[,i], breaks=30, col="magenta") }
}

## set parms for priors and Z
Prior1 = list(ncomp=5)
keep = 5
Mcmc1 = list(R=R, keep=keep)
Data1 = list(p=p, lgtdata=simlgtdata, Z=Z)

## fit model without sign constraints
out1 = rhierMnlRwMixture(Data=Data1, Prior=Prior1, Mcmc=Mcmc1)

cat("Summary of Delta draws", fill=TRUE)
summary(out1$Deltadraw, tvalues=as.vector(Delta))

cat("Summary of Normal Mixture Distribution", fill=TRUE)
summary(out1$nmix)

## plotting examples
if(0) {
  plot(out1$betadraw)
  plot(out1$nmix)
}

## fit model with constraint that beta_i2 < 0 forall i
Prior2 = list(ncomp=5, SignRes=c(0,-1,0))
out2 = rhierMnlRwMixture(Data=Data1, Prior=Prior2, Mcmc=Mcmc1)

cat("Summary of Delta draws", fill=TRUE)
summary(out2$Deltadraw, tvalues=as.vector(Delta))

cat("Summary of Normal Mixture Distribution", fill=TRUE)
summary(out2$nmix)

## plotting examples
if(0) {
  plot(out2$betadraw)
  plot(out2$nmix)
}
}

\keyword{models}
