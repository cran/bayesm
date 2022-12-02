\name{rivGibbs}
\alias{rivGibbs}
\concept{Instrumental Variables}
\concept{Gibbs Sampler}
\concept{bayes}
\concept{endogeneity}
\concept{simultaneity}
\concept{MCMC}

\title{Gibbs Sampler for Linear "IV" Model}

\description{
\code{rivGibbs} is a Gibbs Sampler for a linear structural equation with an arbitrary number of instruments.
}

\usage{rivGibbs(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y, x, w, z)}
  \item{Prior}{list(md, Ad, mbg, Abg, nu, V)}
  \item{Mcmc }{list(R, keep, nprint)}
}

\details{
\subsection{Model and Priors}{
  \eqn{x = z'\delta + e1} \cr
  \eqn{y = \beta*x + w'\gamma + e2} \cr
  \eqn{e1,e2} \eqn{\sim}{~} \eqn{N(0, \Sigma)} 
  
  Note: if intercepts are desired in either equation, include vector of ones in \eqn{z} or \eqn{w}

  \eqn{\delta}            \eqn{\sim}{~} \eqn{N(md, Ad^{-1})} \cr
  \eqn{vec(\beta,\gamma)} \eqn{\sim}{~} \eqn{N(mbg, Abg^{-1})} \cr
  \eqn{\Sigma}            \eqn{\sim}{~} \eqn{IW(nu, V)}
}
\subsection{Argument Details}{
  \emph{\code{Data = list(y, x, w, z)}}
  \tabular{ll}{
    \code{y:      } \tab \eqn{n x 1} vector of obs on LHS variable in structural equation \cr
    \code{x:      } \tab \eqn{n x 1} vector of obs on "endogenous" variable in structural equation \cr
    \code{w:      } \tab \eqn{n x j} matrix of obs on "exogenous" variables in the structural equation \cr
    \code{z:      } \tab \eqn{n x p} matrix of obs on instruments
  }
  \emph{\code{Prior = list(md, Ad, mbg, Abg, nu, V)} [optional]}
  \tabular{ll}{
    \code{md:     } \tab \eqn{p}-length prior mean of delta (def: 0) \cr
    \code{Ad:     } \tab \eqn{p x p} PDS prior precision matrix for prior on delta (def: 0.01*I) \cr
    \code{mbg:    } \tab \eqn{(j+1)}-length prior mean vector for prior on beta,gamma (def: 0) \cr
    \code{Abg:    } \tab \eqn{(j+1)x(j+1)} PDS prior precision matrix for prior on beta,gamma (def: 0.01*I) \cr
    \code{nu:     } \tab d.f. parameter for Inverted Wishart prior on Sigma (def: 5) \cr
    \code{V:      } \tab \eqn{2 x 2} location matrix for Inverted Wishart prior on Sigma (def: nu*I)
  }
  \emph{\code{Mcmc = list(R, keep, nprint)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:      } \tab number of MCMC draws \cr
    \code{keep:   } \tab MCMC thinning parameter: keep every \code{keep}th draw (def: 1) \cr
    \code{nprint: } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print)
  }
}
}

\value{
  A list containing:
  \item{deltadraw }{ \eqn{R/keep x p} matrix of delta draws}
  \item{betadraw  }{ \eqn{R/keep x 1} vector of beta draws}
  \item{gammadraw }{ \eqn{R/keep x j} matrix of gamma draws}
  \item{Sigmadraw }{ \eqn{R/keep x 4} matrix of Sigma draws -- each row is the vector form of Sigma}
}

\author{Rob McCulloch and Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{ For further discussion, see Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

simIV = function(delta, beta, Sigma, n, z, w, gamma) {
  eps = matrix(rnorm(2*n),ncol=2) \%*\% chol(Sigma)
  x = z\%*\%delta + eps[,1]
  y = beta*x + eps[,2] + w\%*\%gamma
  list(x=as.vector(x), y=as.vector(y)) 
  }
  
n = 200
p=1 # number of instruments
z = cbind(rep(1,n), matrix(runif(n*p),ncol=p))
w = matrix(1,n,1)
rho = 0.8
Sigma = matrix(c(1,rho,rho,1), ncol=2)
delta = c(1,4)
beta = 0.5
gamma = c(1)
simiv = simIV(delta, beta, Sigma, n, z, w, gamma)

Data1 = list(); Data1$z = z; Data1$w=w; Data1$x=simiv$x; Data1$y=simiv$y
Mcmc1=list(); Mcmc1$R = R; Mcmc1$keep=1

out = rivGibbs(Data=Data1, Mcmc=Mcmc1)

cat("Summary of Beta draws", fill=TRUE)
summary(out$betadraw, tvalues=beta)

cat("Summary of Sigma draws", fill=TRUE)
summary(out$Sigmadraw, tvalues=as.vector(Sigma[upper.tri(Sigma,diag=TRUE)]))

## plotting examples
if(0){plot(out$betadraw)}
}

\keyword{models}
