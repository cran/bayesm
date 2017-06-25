\name{rnmixGibbs}
\alias{rnmixGibbs}
\concept{bayes}
\concept{MCMC}
\concept{normal mixtures}
\concept{Gibbs Sampling}

\title{Gibbs Sampler for Normal Mixtures}

\description{
\code{rnmixGibbs} implements a Gibbs Sampler for normal mixtures.
}

\usage{rnmixGibbs(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y)}
  \item{Prior}{list(Mubar, A, nu, V, a, ncomp)}
  \item{Mcmc }{list(R, keep, nprint, Loglike)}
}

\details{
\subsection{Model and Priors}{
  \eqn{y_i} \eqn{\sim}{~} \eqn{N(\mu_{ind_i}, \Sigma_{ind_i})} \cr
  ind \eqn{\sim}{~} iid multinomial(p) with \eqn{p} an \eqn{ncomp x 1} vector of probs

  \eqn{\mu_j}    \eqn{\sim}{~} \eqn{N(mubar, \Sigma_j (x) A^{-1})} with \eqn{mubar=vec(Mubar)} \cr
  \eqn{\Sigma_j} \eqn{\sim}{~} \eqn{IW(nu, V)} \cr
  Note: this is the natural conjugate prior -- a special case of multivariate regression \cr
  
  \eqn{p} \eqn{\sim}{~} Dirchlet(a)
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(y)}}
  \tabular{ll}{
    \code{y:       } \tab \eqn{n x k} matrix of data (rows are obs)
    }
  \emph{\code{Prior = list(Mubar, A, nu, V, a, ncomp)} [only \code{ncomp} required]}
  \tabular{ll}{
    \code{Mubar:   } \tab \eqn{1 x k} vector with prior mean of normal component means (def: 0) \cr
    \code{A:       } \tab \eqn{1 x 1} precision parameter for prior on mean of normal component (def: 0.01) \cr
    \code{nu:      } \tab d.f. parameter for prior on Sigma (normal component cov matrix) (def: k+3) \cr
    \code{V:       } \tab \eqn{k x k} location matrix of IW prior on Sigma (def: nu*I) \cr
    \code{a:       } \tab \eqn{ncomp x 1} vector of Dirichlet prior parameters (def: \code{rep(5,ncomp)}) \cr
    \code{ncomp:   } \tab number of normal components to be included
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, Loglike)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:       } \tab number of MCMC draws \cr
    \code{keep:    } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:  } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{LogLike: } \tab logical flag for whether to compute the log-likelihood (def: \code{FALSE})
    }
}
\subsection{\code{nmix} Details}{
  \code{nmix} is a list with 3 components. Several functions in the \code{bayesm} package that involve a Dirichlet Process or mixture-of-normals return \code{nmix}. Across these functions, a common structure is used for \code{nmix} in order to utilize generic summary and plotting functions. 
  \tabular{ll}{
  \code{probdraw:} \tab \eqn{ncomp x R/keep} matrix that reports the probability that each draw came from a particular component \cr
  \code{zdraw:   } \tab \eqn{R/keep x nobs} matrix that indicates which component each draw is assigned to \cr
  \code{compdraw:} \tab A list of \eqn{R/keep} lists of \eqn{ncomp} lists. Each of the inner-most lists has 2 elemens: a vector of draws for \code{mu} and a matrix of draws for the Cholesky root of \code{Sigma}.
  }
}
}

\value{
  A list containing:
  \item{ll   }{ \eqn{R/keep x 1} vector of log-likelihood values}
  \item{nmix }{ a list containing: \code{probdraw}, \code{zdraw}, \code{compdraw} (see \dQuote{\code{nmix} Details} section)}
}

\note{
In this model, the component normal parameters are not-identified due to label-switching. However, the fitted mixture of normals density is identified as it is invariant to label-switching. See chapter 5 of Rossi et al below for details. 
  
Use \code{eMixMargDen} or \code{momMix} to compute posterior expectation or distribution of various identified parameters.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. \cr \url{http://www.perossi.org/home/bsm-1}}

\seealso{ \code{\link{rmixture}}, \code{\link{rmixGibbs}} ,\code{\link{eMixMargDen}}, \code{\link{momMix}},
 \code{\link{mixDen}}, \code{\link{mixDenBi}}}

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

dim = 5
k = 3   # dimension of simulated data and number of "true" components
sigma = matrix(rep(0.5,dim^2), nrow=dim)
diag(sigma) = 1
sigfac = c(1,1,1)
mufac = c(1,2,3)
compsmv = list()
for(i in 1:k) compsmv[[i]] = list(mu=mufac[i]*1:dim, sigma=sigfac[i]*sigma)

# change to "rooti" scale
comps = list() 
for(i in 1:k) comps[[i]] = list(mu=compsmv[[i]][[1]], rooti=solve(chol(compsmv[[i]][[2]])))
pvec = (1:k) / sum(1:k)

nobs = 500
dm = rmixture(nobs, pvec, comps)

Data1 = list(y=dm$x)
ncomp = 9
Prior1 = list(ncomp=ncomp)
Mcmc1 = list(R=R, keep=1)

out = rnmixGibbs(Data=Data1, Prior=Prior1, Mcmc=Mcmc1)

cat("Summary of Normal Mixture Distribution", fill=TRUE)
summary(out$nmix)

tmom = momMix(matrix(pvec,nrow=1), list(comps))
mat = rbind(tmom$mu, tmom$sd)
cat(" True Mean/Std Dev", fill=TRUE)
print(mat)

## plotting examples
if(0){plot(out$nmix,Data=dm$x)}
}

\keyword{multivariate}
