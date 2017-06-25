\name{rivDP}
\alias{rivDP}
\concept{Instrumental Variables}
\concept{Gibbs Sampler}
\concept{Dirichlet Process}
\concept{bayes}
\concept{endogeneity}
\concept{simultaneity}
\concept{MCMC}

\title{Linear "IV" Model with DP Process Prior for Errors}

\description{
\code{rivDP} is a Gibbs Sampler for a linear structural equation with an arbitrary number of instruments. \code{rivDP} uses a mixture-of-normals for the structural and reduced form equations implemented with a Dirichlet Process prior.
}

\usage{rivDP(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y, x, w, z)}
  \item{Prior}{list(md, Ad, mbg, Abg, lambda, Prioralpha, lambda_hyper)}
  \item{Mcmc }{list(R, keep, nprint, maxuniq, SCALE, gridsize)}
}

\details{
\subsection{Model and Priors}{
  \eqn{x = z'\delta + e1} \cr
  \eqn{y = \beta*x + w'\gamma + e2} \cr
  \eqn{e1,e2} \eqn{\sim}{~} \eqn{N(\theta_{i})} where \eqn{\theta_{i}} represents \eqn{\mu_{i}, \Sigma_{i}}
  
  Note: Error terms have non-zero means. 
  DO NOT include intercepts in the \eqn{z} or \eqn{w} matrices.  
  This is different from \code{rivGibbs} which requires intercepts to be included explicitly.

  \eqn{\delta} \eqn{\sim}{~} \eqn{N(md, Ad^{-1})} \cr
  \eqn{vec(\beta, \gamma)} \eqn{\sim}{~} \eqn{N(mbg, Abg^{-1})} \cr
  \eqn{\theta_{i}} \eqn{\sim}{~} \eqn{G} \cr
  \eqn{G} \eqn{\sim}{~} \eqn{DP(alpha, G_0)} 
 
  \eqn{alpha} \eqn{\sim}{~} \eqn{(1-(alpha-alpha_{min})/(alpha_{max}-alpha{min}))^{power}} \cr 
  where \eqn{alpha_{min}} and \eqn{alpha_{max}} are set using the arguments in the reference below.  
  It is highly recommended that you use the default values for the hyperparameters of the prior on alpha.
 
  \eqn{G_0} is the natural conjugate prior for \eqn{(\mu,\Sigma)}: 
  \eqn{\Sigma} \eqn{\sim}{~} \eqn{IW(nu, vI)} and  \eqn{\mu|\Sigma} \eqn{\sim}{~} \eqn{N(0, \Sigma(x) a^{-1})} \cr
  These parameters are collected together in the list \eqn{\lambda}.  
  It is highly recommended that you use the default settings for these hyper-parameters.\cr
  
  \eqn{\lambda(a, nu, v):}\cr
    \eqn{a}  \eqn{\sim}{~} uniform[alim[1], alimb[2]]\cr
    \eqn{nu} \eqn{\sim}{~} dim(data)-1 + exp(z) \cr
    \eqn{z}  \eqn{\sim}{~} uniform[dim(data)-1+nulim[1], nulim[2]]\cr
    \eqn{v}  \eqn{\sim}{~} uniform[vlim[1], vlim[2]]
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(y, x, w, z)}}
  \tabular{ll}{
    \code{y: } \tab \eqn{n x 1} vector of obs on LHS variable in structural equation \cr
    \code{x: } \tab \eqn{n x 1} vector of obs on "endogenous" variable in structural equation \cr
    \code{w: } \tab \eqn{n x j} matrix of obs on "exogenous" variables in the structural equation \cr
    \code{z: } \tab \eqn{n x p} matrix of obs on instruments
  }
  \emph{\code{Prior = list(md, Ad, mbg, Abg, lambda, Prioralpha, lambda_hyper)} [optional]}
  \tabular{ll}{
    \code{md:  } \tab \eqn{p}-length prior mean of delta (def: 0) \cr
    \code{Ad:  } \tab \eqn{p x p} PDS prior precision matrix for prior on delta (def: 0.01*I) \cr
    \code{mbg: } \tab \eqn{(j+1)}-length prior mean vector for prior on beta,gamma (def: 0) \cr
    \code{Abg: } \tab \eqn{(j+1)x(j+1)} PDS prior precision matrix for prior on beta,gamma (def: 0.01*I) \cr
    \code{Prioralpha:} \tab \code{list(Istarmin, Istarmax, power)} \cr
    \code{$Istarmin: } \tab  is expected number of components at lower bound of support of alpha (def: 1) \cr
    \code{$Istarmax: } \tab  is expected number of components at upper bound of support of alpha (def: \code{floor(0.1*length(y))}) \cr
    \code{$power:    } \tab  is the power parameter for alpha prior (def: 0.8)  \cr
    \code{lambda_hyper:} \tab \code{list(alim, nulim, vlim)} \cr
    \code{$alim:     } \tab  defines support of a distribution (def: \code{c(0.01, 10)}) \cr
    \code{$nulim:    } \tab  defines support of nu distribution (def: \code{c(0.01, 3)}) \cr
    \code{$vlim:     } \tab  defines support of v distribution (def: \code{c(0.1, 4)}) 
  }
  \emph{\code{Mcmc  = list(R, keep, nprint, maxuniq, SCALE, gridsize)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:        } \tab number of MCMC draws \cr
    \code{keep:     } \tab MCMC thinning parameter: keep every keepth draw (def: 1) \cr
    \code{nprint:   } \tab print the estimated time remaining for every nprint'th draw (def: 100, set to 0 for no print) \cr
    \code{maxuniq:  } \tab storage constraint on the number of unique components (def: 200) \cr
    \code{SCALE:    } \tab scale data (def: \code{TRUE}) \cr
    \code{gridsize: } \tab gridsize parameter for alpha draws (def: 20)
  }
}
\subsection{\code{nmix} Details}{
  \code{nmix} is a list with 3 components. Several functions in the \code{bayesm} package that involve a Dirichlet Process or mixture-of-normals return \code{nmix}. Across these functions, a common structure is used for \code{nmix} in order to utilize generic summary and plotting functions. 
  \tabular{ll}{
  \code{probdraw:} \tab \eqn{ncomp x R/keep} matrix that reports the probability that each draw came from a particular component (here, a one-column matrix of 1s) \cr
  \code{zdraw:   } \tab \eqn{R/keep x nobs} matrix that indicates which component each draw is assigned to (here, null) \cr
  \code{compdraw:} \tab A list of \eqn{R/keep} lists of \eqn{ncomp} lists. Each of the inner-most lists has 2 elemens: a vector of draws for \code{mu} and a matrix of draws for the Cholesky root of \code{Sigma}.
  }
}
}

\value{
  A list containing:
  \item{\code{deltadraw }}{ \eqn{R/keep x p} array of delta draws}
  \item{\code{betadraw  }}{ \eqn{R/keep x 1} vector of beta draws}
  \item{\code{alphadraw }}{ \eqn{R/keep x 1} vector of draws of Dirichlet Process tightness parameter}
  \item{\code{Istardraw }}{ \eqn{R/keep x 1} vector of draws of the number of unique normal components}
  \item{\code{gammadraw }}{ \eqn{R/keep x j} array of gamma draws}
  \item{\code{nmix      }}{ a list containing: \code{probdraw}, \code{zdraw}, \code{compdraw} (see \dQuote{\code{nmix} Details} section)}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{
For further discussion, see "A Semi-Parametric Bayesian Approach to the Instrumental Variable Problem," by Conley, Hansen, McCulloch and Rossi, \emph{Journal of Econometrics} (2008). 

See also, Chapter 4, \emph{Bayesian Non- and Semi-parametric Methods and Applications} by Peter Rossi.
}

\seealso{\code{rivGibbs}}

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

## simulate scaled log-normal errors and run
k = 10
delta = 1.5
Sigma = matrix(c(1, 0.6, 0.6, 1), ncol=2)
N = 1000
tbeta = 4
scalefactor = 0.6
root = chol(scalefactor*Sigma)
mu = c(1,1)

## compute interquartile ranges
ninterq = qnorm(0.75) - qnorm(0.25)
error = matrix(rnorm(100000*2), ncol=2)\%*\%root
error = t(t(error)+mu)
Err = t(t(exp(error))-exp(mu+0.5*scalefactor*diag(Sigma)))
lnNinterq = quantile(Err[,1], prob=0.75) - quantile(Err[,1], prob=0.25)

## simulate data
error = matrix(rnorm(N*2), ncol=2)\%*\%root
error = t(t(error)+mu)
Err = t(t(exp(error))-exp(mu+0.5*scalefactor*diag(Sigma)))

## scale appropriately  
Err[,1] = Err[,1]*ninterq/lnNinterq
Err[,2] = Err[,2]*ninterq/lnNinterq
z = matrix(runif(k*N), ncol=k)
x = z\%*\%(delta*c(rep(1,k))) + Err[,1]
y = x*tbeta + Err[,2]

## specify data input and mcmc parameters
Data = list(); 
Data$z = z
Data$x = x
Data$y = y

Mcmc = list()
Mcmc$maxuniq = 100
Mcmc$R = R
end = Mcmc$R

out = rivDP(Data=Data, Mcmc=Mcmc)

cat("Summary of Beta draws", fill=TRUE)
summary(out$betadraw, tvalues=tbeta)

## plotting examples
if(0){
  plot(out$betadraw, tvalues=tbeta)
  plot(out$nmix)  # plot "fitted" density of the errors
}
}

\keyword{models}
