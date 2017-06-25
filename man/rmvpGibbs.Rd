\name{rmvpGibbs}
\alias{rmvpGibbs}

\concept{bayes}
\concept{multivariate probit}
\concept{MCMC}
\concept{Gibbs Sampling}

\title{Gibbs Sampler for Multivariate Probit}

\description{
\code{rmvpGibbs} implements the Edwards/Allenby Gibbs Sampler for the multivariate probit model.
}
  
\usage{rmvpGibbs(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y, X, p)}
  \item{Prior}{list(betabar, A, nu, V)}
  \item{Mcmc }{list(R, keep, nprint, beta0 ,sigma0)}
}

\details{
\subsection{Model and Priors}{
  \eqn{w_i = X_i\beta + e} with \eqn{e} \eqn{\sim}{~} N(0,\eqn{\Sigma}). Note: \eqn{w_i} is \eqn{p x 1}. \cr
  \eqn{y_{ij} = 1} if \eqn{w_{ij} > 0}, else \eqn{y_i = 0}.  \eqn{j = 1, \ldots, p} \cr
  
  beta and Sigma are not identifed. Correlation matrix and the betas divided by the appropriate standard deviation are. 
  See reference or example below for details.

  \eqn{\beta}  \eqn{\sim}{~} \eqn{N(betabar, A^{-1})} \cr
  \eqn{\Sigma} \eqn{\sim}{~} \eqn{IW(nu, V)} \cr

  To make \eqn{X} matrix use \code{createX}
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(y, X, p)}}
  \tabular{ll}{
    \code{X:       } \tab \eqn{n*p x k} Design Matrix \cr
    \code{y:       } \tab \eqn{n*p x 1} vector of 0/1 outcomes \cr
    \code{p:       } \tab dimension of multivariate probit
    }
  \emph{\code{Prior = list(betabar, A, nu, V)} [optional]}
  \tabular{ll}{
    \code{betabar: } \tab \eqn{k x 1} prior mean (def: 0) \cr
    \code{A:       } \tab \eqn{k x k} prior precision matrix (def: 0.01*I) \cr
    \code{nu:      } \tab d.f. parameter for Inverted Wishart prior (def: (p-1)+3) \cr
    \code{V:       } \tab PDS location parameter for Inverted Wishart prior (def: nu*I) 
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, beta0 ,sigma0)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:       } \tab number of MCMC draws \cr
    \code{keep:    } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:  } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{beta0:   } \tab initial value for beta \cr
    \code{sigma0:  } \tab initial value for sigma
    }
}
}

\value{
  A list containing: 
  \item{betadraw  }{\eqn{R/keep x k} matrix of betadraws}
  \item{sigmadraw }{\eqn{R/keep x p*p} matrix of sigma draws -- each row is the vector form of sigma}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 4, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. \cr \url{http://www.perossi.org/home/bsm-1}}

\seealso{ \code{\link{rmnpGibbs}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

simmvp = function(X, p, n, beta, sigma) {
  w = as.vector(crossprod(chol(sigma),matrix(rnorm(p*n),ncol=n))) + X\%*\%beta
  y = ifelse(w<0, 0, 1)
  return(list(y=y, X=X, beta=beta, sigma=sigma))
}

p = 3
n = 500
beta = c(-2,0,2)
Sigma = matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), ncol=3)
k = length(beta)
I2 = diag(rep(1,p))
xadd = rbind(I2)
for(i in 2:n) { xadd=rbind(xadd,I2) }
X = xadd

simout = simmvp(X,p,500,beta,Sigma)

Data1 = list(p=p, y=simout$y, X=simout$X)
Mcmc1 = list(R=R, keep=1)

out = rmvpGibbs(Data=Data1, Mcmc=Mcmc1)

ind = seq(from=0, by=p, length=k)
inda = 1:3
ind = ind + inda
cat(" Betadraws ", fill=TRUE)
betatilde = out$betadraw / sqrt(out$sigmadraw[,ind])
attributes(betatilde)$class = "bayesm.mat"
summary(betatilde, tvalues=beta/sqrt(diag(Sigma)))

rdraw = matrix(double((R)*p*p), ncol=p*p)
rdraw = t(apply(out$sigmadraw, 1, nmat))
attributes(rdraw)$class = "bayesm.var"
tvalue = nmat(as.vector(Sigma))
dim(tvalue) = c(p,p)
tvalue = as.vector(tvalue[upper.tri(tvalue,diag=TRUE)])
cat(" Draws of Correlation Matrix ", fill=TRUE)
summary(rdraw, tvalues=tvalue)

## plotting examples
if(0){plot(betatilde, tvalues=beta/sqrt(diag(Sigma)))}
}

\keyword{models}
\keyword{multivariate}
