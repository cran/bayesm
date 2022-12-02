\name{rmnlIndepMetrop}
\alias{rmnlIndepMetrop}
\concept{MCMC}
\concept{multinomial logit}
\concept{Metropolis algorithm}
\concept{bayes}

\title{MCMC Algorithm for Multinomial Logit Model}

\description{
\code{rmnlIndepMetrop} implements Independence Metropolis algorithm for the multinomial logit (MNL) model.
}

\usage{rmnlIndepMetrop(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y, X, p)}
  \item{Prior}{list(A, betabar)}
  \item{Mcmc }{list(R, keep, nprint, nu)}
}

\details{
\subsection{Model and Priors}{
  y \eqn{\sim}{~} MNL(X, \eqn{\beta}) \cr
  \eqn{Pr(y=j) = exp(x_j'\beta)/\sum_k{e^{x_k'\beta}}}

  \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar, A^{-1})}
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(y, X, p)}}
  \tabular{ll}{
    \code{y:       } \tab \eqn{n x 1} vector of multinomial outcomes (1, \ldots, p) \cr
    \code{X:       } \tab \eqn{n*p x k} matrix \cr
    \code{p:       } \tab number of alternatives
  }
  \emph{\code{Prior = list(A, betabar)} [optional]}
  \tabular{ll}{
    \code{A:       } \tab \eqn{k x k} prior precision matrix (def: 0.01*I) \cr
    \code{betabar: } \tab \eqn{k x 1} prior mean (def: 0)
  }
  \emph{\code{Mcmc = list(R, keep, nprint, nu)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:       } \tab number of MCMC draws \cr
    \code{keep:    } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:  } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{nu:      } \tab d.f. parameter for independent t density (def: 6)
  }
}
}

\value{
  A list containing: 
  \item{betadraw }{\eqn{R/keep x k} matrix of beta draws}
  \item{loglike  }{\eqn{R/keep x 1} vector of log-likelihood values evaluated at each draw}
  \item{acceptr  }{acceptance rate of Metropolis draws}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{rhierMnlRwMixture}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

simmnl = function(p, n, beta) {
  ## note: create X array with 2 alt.spec vars
    k = length(beta)
    X1 = matrix(runif(n*p,min=-1,max=1), ncol=p)
    X2 = matrix(runif(n*p,min=-1,max=1), ncol=p)
    X = createX(p, na=2, nd=NULL, Xd=NULL, Xa=cbind(X1,X2), base=1)
    Xbeta = X\%*\%beta 
  ## now do probs
    p = nrow(Xbeta) / n
    Xbeta = matrix(Xbeta, byrow=TRUE, ncol=p)
    Prob = exp(Xbeta)
    iota = c(rep(1,p))
    denom = Prob\%*\%iota
    Prob = Prob / as.vector(denom)
  ## draw y
    y = vector("double",n)
    ind = 1:p
    for (i in 1:n) { 
      yvec = rmultinom(1, 1, Prob[i,])
      y[i] = ind\%*\%yvec 
    }
  return(list(y=y, X=X, beta=beta, prob=Prob))
}

n = 200
p = 3
beta = c(1, -1, 1.5, 0.5)

simout = simmnl(p,n,beta)

Data1 = list(y=simout$y, X=simout$X, p=p)
Mcmc1 = list(R=R, keep=1)

out = rmnlIndepMetrop(Data=Data1, Mcmc=Mcmc1)

cat("Summary of beta draws", fill=TRUE)
summary(out$betadraw, tvalues=beta)

## plotting examples
if(0){plot(out$betadraw)}
}

\keyword{models}
