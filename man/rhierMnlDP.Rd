\name{rhierMnlDP}
\alias{rhierMnlDP}
\concept{bayes}
\concept{MCMC}
\concept{Multinomial Logit}
\concept{normal mixture}
\concept{Dirichlet Process Prior}
\concept{heterogeneity}
\concept{hierarchical models}

\title{MCMC Algorithm for Hierarchical Multinomial Logit with Dirichlet Process Prior Heterogeneity}

\description{
\code{rhierMnlDP} is a MCMC algorithm for a hierarchical multinomial logit with a Dirichlet Process prior for the distribution of heteorogeneity.  A base normal model is used so that the DP can be interpreted as allowing for a mixture of normals with as many components as there are panel units. This is a hybrid Gibbs Sampler with a RW Metropolis step for the MNL coefficients for each panel unit.  This procedure can be interpreted as a Bayesian semi-parameteric method in the sense that the DP prior can accomodate heterogeniety of an unknown form.
}

\usage{rhierMnlDP(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(lgtdata, Z, p)}
  \item{Prior}{list(deltabar, Ad, Prioralpha, lambda_hyper)}
  \item{Mcmc }{list(R, keep, nprint, s, w, maxunique, gridsize)}
}

\details{
\subsection{Model and Priors}{
  \eqn{y_i} \eqn{\sim}{~} \eqn{MNL(X_i, \beta_i)} with \eqn{i = 1, \ldots, length(lgtdata)} and where \eqn{\theta_i} is \eqn{nvar x 1}

  \eqn{\beta_i = Z\Delta}[i,] + \eqn{u_i} \cr
  Note:  Z\eqn{\Delta} is the matrix \eqn{Z * \Delta}; [i,] refers to \eqn{i}th row of this product \cr
  Delta is an \eqn{nz x nvar} matrix 

  \eqn{\beta_i} \eqn{\sim}{~} \eqn{N(\mu_i, \Sigma_i)}

  \eqn{\theta_i = (\mu_i, \Sigma_i)} \eqn{\sim}{~} \eqn{DP(G_0(\lambda), alpha)}\cr
  
  \eqn{G_0(\lambda):}\cr
  \eqn{\mu_i | \Sigma_i} \eqn{\sim}{~} \eqn{N(0, \Sigma_i (x) a^{-1})}\cr
  \eqn{\Sigma_i} \eqn{\sim}{~} \eqn{IW(nu, nu*v*I)}\cr
  \eqn{delta = vec(\Delta)} \eqn{\sim}{~} \eqn{N(deltabar, A_d^{-1})}\cr
  
  \eqn{\lambda(a, nu, v):}\cr
  \eqn{a} \eqn{\sim}{~} uniform[alim[1], alimb[2]]\cr
  \eqn{nu} \eqn{\sim}{~}  dim(data)-1 + exp(z) \cr
  \eqn{z} \eqn{\sim}{~}  uniform[dim(data)-1+nulim[1], nulim[2]]\cr
  \eqn{v} \eqn{\sim}{~} uniform[vlim[1], vlim[2]]
  
  \eqn{alpha} \eqn{\sim}{~} \eqn{(1-(alpha-alphamin) / (alphamax-alphamin))^{power}} \cr
  alpha = alphamin then expected number of components = \code{Istarmin} \cr
  alpha = alphamax then expected number of components = \code{Istarmax}
  
  \eqn{Z} should NOT include an intercept and is centered for ease of interpretation. The mean of each of the \code{nlgt} \eqn{\beta}s is the mean of the normal mixture.  Use \code{summary()} to compute this mean from the \code{compdraw} output.
  
  We parameterize the prior on \eqn{\Sigma_i} such that \eqn{mode(\Sigma) = nu/(nu+2) vI}. The support of nu enforces a non-degenerate IW density; \eqn{nulim[1] > 0}.

  The default choices of alim, nulim, and vlim determine the location and approximate size of candidate "atoms" or possible normal components. The defaults are sensible given a reasonable scaling of the X variables. You want to insure that alim is set for a wide enough range of values (remember a is a precision parameter) and the v is big enough to propose Sigma matrices wide enough to cover the data range.  

  A careful analyst should look at the posterior distribution of a, nu, v to make sure that the support is set correctly in alim, nulim, vlim.  In other words, if we see the posterior bunched up at one end of these support ranges, we should widen the range and rerun.  

  If you want to force the procedure to use many small atoms, then set nulim to consider only large values and set vlim to consider only small scaling constants.  Set alphamax to a large number.  This will create a very "lumpy" density estimate somewhat like the classical Kernel density estimates. Of course, this is not advised if you have a prior belief that densities are relatively smooth.
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(lgtdata, Z, p)} [\code{Z} optional]}
  \tabular{ll}{
    \code{lgtdata:        } \tab list of lists with each cross-section unit MNL data \cr
    \code{lgtdata[[i]]$y: } \tab \eqn{n_i x 1} vector of multinomial outcomes (1, \ldots, m) \cr
    \code{lgtdata[[i]]$X: } \tab \eqn{n_i x nvar} design matrix for \eqn{i}th unit \cr
    \code{Z:              } \tab \eqn{nreg x nz} matrix of unit characteristics (def: vector of ones) \cr
    \code{p:              } \tab number of choice alternatives
    }
  \emph{\code{Prior = list(deltabar, Ad, Prioralpha, lambda_hyper)} [optional]}
  \tabular{ll}{
    \code{deltabar:       } \tab \eqn{nz*nvar x 1} vector of prior means (def: 0) \cr
    \code{Ad:             } \tab prior precision matrix for vec(D) (def: 0.01*I) \cr
    \code{Prioralpha:     } \tab \code{list(Istarmin, Istarmax, power)} \cr
      \code{$Istarmin:    } \tab expected number of components at lower bound of support of alpha def(1) \cr
      \code{$Istarmax:    } \tab expected number of components at upper bound of support of alpha (def: min(50, 0.1*nlgt)) \cr
      \code{$power:       } \tab power parameter for alpha prior (def: 0.8)  \cr
    \code{lambda_hyper:   } \tab \code{list(alim, nulim, vlim)}  \cr
      \code{$alim:        } \tab defines support of a distribution (def: \code{c(0.01, 2)}) \cr
      \code{$nulim:       } \tab defines support of nu distribution (def: \code{c(0.001, 3)}) \cr
      \code{$vlim:        } \tab defines support of v distribution (def: \code{c(0.1, 4)}) 
      }
  \emph{\code{Mcmc  = list(R, keep, nprint, s, w, maxunique, gridsize)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:              } \tab number of MCMC draws \cr
    \code{keep:           } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:         } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{s:              } \tab scaling parameter for RW Metropolis (def: 2.93/\code{sqrt(nvar)}) \cr
    \code{w:              } \tab fractional likelihood weighting parameter (def: 0.1) \cr
    \code{maxuniq:        } \tab storage constraint on the number of unique components (def: 200) \cr
    \code{gridsize:       } \tab number of discrete points for hyperparameter priors, (def: 20)
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
  \item{Deltadraw }{ \eqn{R/keep x nz*nvar} matrix of draws of Delta, first row is initial value}
  \item{betadraw  }{ \eqn{nlgt x nvar x R/keep} array of draws of betas}
  \item{nmix      }{ a list containing: \code{probdraw}, \code{zdraw}, \code{compdraw} (see \dQuote{\code{nmix} Details} section)}
  \item{adraw     }{ \eqn{R/keep} draws of hyperparm a}
  \item{vdraw     }{ \eqn{R/keep} draws of hyperparm v}
  \item{nudraw    }{ \eqn{R/keep} draws of hyperparm nu}
  \item{Istardraw }{ \eqn{R/keep} draws of number of unique components}
  \item{alphadraw }{ \eqn{R/keep} draws of number of DP tightness parameter}
  \item{loglike   }{ \eqn{R/keep} draws of log-likelihood}
}

\note{
As is well known, Bayesian density estimation involves computing the predictive distribution of a "new" unit parameter, \eqn{\theta_{n+1}} (here "n"=nlgt). This is done by averaging the normal base distribution over draws from the distribution of \eqn{\theta_{n+1}} given \eqn{\theta_1}, ..., \eqn{\theta_n}, alpha, lambda, data. To facilitate this, we store those draws from the predictive distribution of \eqn{\theta_{n+1}} in a list structure compatible with other \code{bayesm} routines that implement a finite mixture of normals.
  
Large \code{R} values may be required (>20,000).
} 

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. \cr \url{http://www.perossi.org/home/bsm-1}}
  
\seealso{ \code{\link{rhierMnlRwMixture}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=20000} else {R=10}
set.seed(66)

p = 3                                # num of choice alterns
ncoef = 3  
nlgt = 300                           # num of cross sectional units
nz = 2
Z = matrix(runif(nz*nlgt),ncol=nz)
Z = t(t(Z)-apply(Z,2,mean))          # demean Z
ncomp = 3                            # no of mixture components
Delta=matrix(c(1,0,1,0,1,2),ncol=2)

comps = NULL
comps[[1]] = list(mu=c(0,-1,-2),   rooti=diag(rep(2,3)))
comps[[2]] = list(mu=c(0,-1,-2)*2, rooti=diag(rep(2,3)))
comps[[3]] = list(mu=c(0,-1,-2)*4, rooti=diag(rep(2,3)))
pvec=c(0.4, 0.2, 0.4)

##  simulate from MNL model conditional on X matrix
simmnlwX = function(n,X,beta) {
  k = length(beta)
  Xbeta = X\%*\%beta
  j = nrow(Xbeta) / n
  Xbeta = matrix(Xbeta, byrow=TRUE, ncol=j)
  Prob = exp(Xbeta)
  iota = c(rep(1,j))
  denom = Prob\%*\%iota
  Prob = Prob / as.vector(denom)
  y = vector("double", n)
  ind = 1:j
  for (i in 1:n) {
  yvec = rmultinom(1, 1, Prob[i,])
  y[i] = ind\%*\%yvec}
  return(list(y=y, X=X, beta=beta, prob=Prob))
}

## simulate data with a mixture of 3 normals
simlgtdata = NULL
ni = rep(50,300)
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
  for(i in 1:nlgt) { bmat[i,] = simlgtdata[[i]]$beta }
  par(mfrow = c(ncoef,1))
  for(i in 1:ncoef) { hist(bmat[,i], breaks=30, col="magenta") }
}

## set Data and Mcmc lists
keep = 5
Mcmc1 = list(R=R, keep=keep)
Data1 = list(p=p, lgtdata=simlgtdata, Z=Z)

out = rhierMnlDP(Data=Data1, Mcmc=Mcmc1)

cat("Summary of Delta draws", fill=TRUE)
summary(out$Deltadraw, tvalues=as.vector(Delta))

## plotting examples
if(0) {
  plot(out$betadraw)
  plot(out$nmix)
}
}

\keyword{models}
