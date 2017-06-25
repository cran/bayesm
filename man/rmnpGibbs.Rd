\name{rmnpGibbs}
\alias{rmnpGibbs}
\concept{bayes}
\concept{multinomial probit}
\concept{MCMC}
\concept{Gibbs Sampling}

\title{Gibbs Sampler for Multinomial Probit}

\description{
\code{rmnpGibbs} implements the McCulloch/Rossi Gibbs Sampler for the multinomial probit model.
}

\usage{rmnpGibbs(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y, X, p)}
  \item{Prior}{list(betabar, A, nu, V)}
  \item{Mcmc }{list(R, keep, nprint, beta0, sigma0)}
}

\details{
\subsection{Model and Priors}{
  \eqn{w_i = X_i\beta + e} with \eqn{e} \eqn{\sim}{~} \eqn{N(0, \Sigma)}. 
  Note: \eqn{w_i} and \eqn{e} are \eqn{(p-1) x 1}.\cr
  \eqn{y_i = j} if \eqn{w_{ij} > max(0,w_{i,-j})} for \eqn{j=1, \ldots, p-1} and 
  where \eqn{w_{i,-j}} means elements of \eqn{w_i} other than the \eqn{j}th. \cr
  \eqn{y_i = p},  if all \eqn{w_i < 0}
    
  \eqn{\beta} is not identified. However, \eqn{\beta}/sqrt(\eqn{\sigma_{11}}) and 
  \eqn{\Sigma}/\eqn{\sigma_{11}} are identified.  See reference or example below for details.

  \eqn{\beta}  \eqn{\sim}{~} \eqn{N(betabar,A^{-1})} \cr
  \eqn{\Sigma} \eqn{\sim}{~} \eqn{IW(nu,V)} \cr
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(y, X, p)}}
  \tabular{ll}{
    \code{y:       } \tab \eqn{n x 1} vector of multinomial outcomes (1, \ldots, p) \cr
    \code{X:       } \tab \eqn{n*(p-1) x k} design matrix. To make \eqn{X} matrix use \code{\link{createX}} with \code{DIFF=TRUE} \cr
    \code{p:       } \tab number of alternatives
    }
  \emph{\code{Prior = list(betabar, A, nu, V)} [optional]}
  \tabular{ll}{
    \code{betabar: } \tab \eqn{k x 1} prior mean (def: 0) \cr
    \code{A:       } \tab \eqn{k x k} prior precision matrix (def: 0.01*I) \cr
    \code{nu:      } \tab d.f. parameter for Inverted Wishart prior (def: (p-1)+3) \cr
    \code{V:       } \tab PDS location parameter for Inverted Wishart prior (def: nu*I) 
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, beta0, sigma0)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:       } \tab number of MCMC draws \cr
    \code{keep:    } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:  } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{beta0:   } \tab initial value for beta (def: 0) \cr
    \code{sigma0:  } \tab initial value for sigma (def: I)
    }
}
}

\value{
  A list containing: 
  \item{betadraw  }{\eqn{R/keep x k} matrix of betadraws}
  \item{sigmadraw }{\eqn{R/keep x (p-1)*(p-1)} matrix of sigma draws -- each row is the vector form of sigma}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 4, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. \cr \url{http://www.perossi.org/home/bsm-1}}

\seealso{ \code{\link{rmvpGibbs}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

simmnp = function(X, p, n, beta, sigma) {
  indmax = function(x) {which(max(x)==x)}
  Xbeta = X\%*\%beta
  w = as.vector(crossprod(chol(sigma),matrix(rnorm((p-1)*n),ncol=n))) + Xbeta
  w = matrix(w, ncol=(p-1), byrow=TRUE)
  maxw = apply(w, 1, max)
  y = apply(w, 1, indmax)
  y = ifelse(maxw < 0, p, y)
  return(list(y=y, X=X, beta=beta, sigma=sigma))
}

p = 3
n = 500
beta = c(-1,1,1,2)
Sigma = matrix(c(1, 0.5, 0.5, 1), ncol=2)
k = length(beta)
X1 = matrix(runif(n*p,min=0,max=2),ncol=p)
X2 = matrix(runif(n*p,min=0,max=2),ncol=p)
X = createX(p, na=2, nd=NULL, Xa=cbind(X1,X2), Xd=NULL, DIFF=TRUE, base=p)

simout = simmnp(X,p,500,beta,Sigma)

Data1 = list(p=p, y=simout$y, X=simout$X)
Mcmc1 = list(R=R, keep=1)

out = rmnpGibbs(Data=Data1, Mcmc=Mcmc1)

cat(" Summary of Betadraws ", fill=TRUE)
betatilde = out$betadraw / sqrt(out$sigmadraw[,1])
attributes(betatilde)$class = "bayesm.mat"
summary(betatilde, tvalues=beta)

cat(" Summary of Sigmadraws ", fill=TRUE)
sigmadraw = out$sigmadraw / out$sigmadraw[,1]
attributes(sigmadraw)$class = "bayesm.var"
summary(sigmadraw, tvalues=as.vector(Sigma[upper.tri(Sigma,diag=TRUE)]))

## plotting examples
if(0){plot(betatilde,tvalues=beta)}
}

\keyword{models}
