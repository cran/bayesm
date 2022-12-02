\name{rordprobitGibbs}
\alias{rordprobitGibbs}
\concept{bayes}
\concept{MCMC}
\concept{probit}
\concept{Gibbs Sampling}

\title{Gibbs Sampler for Ordered Probit}

\description{
\code{rordprobitGibbs} implements a Gibbs Sampler for the ordered probit model with a RW Metropolis step for the cut-offs.
}

\usage{rordprobitGibbs(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(y, X, k)}
  \item{Prior}{list(betabar, A, dstarbar, Ad)}
  \item{Mcmc }{list(R, keep, nprint, s)}
}

\details{
\subsection{Model and Priors}{
  \eqn{z = X\beta + e} with \eqn{e} \eqn{\sim}{~} \eqn{N(0, I)}\cr 
  \eqn{y = k} if c[k] \eqn{\le z <} c[k+1] with \eqn{k = 1,\ldots,K} \cr
  cutoffs = \{c[1], \eqn{\ldots}, c[k+1]\} 

  \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar, A^{-1})} \cr
  \eqn{dstar} \eqn{\sim}{~} \eqn{N(dstarbar, Ad^{-1})}
  
  Be careful in assessing prior parameter \code{Ad}:  0.1 is too small for many applications.
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(y, X, k)}}
  \tabular{ll}{
    \code{y:        } \tab \eqn{n x 1} vector of observations, (\eqn{1, \ldots, k}) \cr
    \code{X:        } \tab \eqn{n x p} Design Matrix \cr
    \code{k:        } \tab the largest possible value of y
    }
  \emph{\code{Prior = list(betabar, A, dstarbar, Ad)} [optional]}
  \tabular{ll}{
    \code{betabar:  } \tab \eqn{p x 1} prior mean (def: 0) \cr
    \code{A:        } \tab \eqn{p x p} prior precision matrix (def: 0.01*I) \cr
    \code{dstarbar: } \tab \eqn{ndstar x 1} prior mean,  where \eqn{ndstar=k-2} (def: 0) \cr
    \code{Ad:       } \tab \eqn{ndstar x ndstar} prior precision matrix (def: I)
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, s)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:        } \tab number of MCMC draws \cr
    \code{keep:     } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:   } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{s:        } \tab scaling parameter for RW Metropolis (def: 2.93/\code{sqrt(p)})
    }
}
}

\value{
  A list containing: 
  \item{betadraw  }{ \eqn{R/keep x p} matrix of betadraws}
  \item{cutdraw   }{ \eqn{R/keep x (k-1)} matrix of cutdraws}
  \item{dstardraw }{ \eqn{R/keep x (k-2)} matrix of dstardraws}
  \item{accept    }{ acceptance rate of Metropolis draws for cut-offs}
}

\note{ 
   set c[1] = -100 and c[k+1] = 100. c[2] is set to 0 for identification. \cr

   The relationship between cut-offs and dstar is: \cr
   c[3] = exp(dstar[1]),  \cr
   c[4] = c[3] + exp(dstar[2]), ...,  \cr
   c[k] = c[k-1] + exp(dstar[k-2]) 
} 

\references{\emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\author{ Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\seealso{ \code{\link{rbprobitGibbs}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

## simulate data for ordered probit model

simordprobit=function(X, betas, cutoff){
  z = X\%*\%betas + rnorm(nobs)   
  y = cut(z, br = cutoff, right=TRUE, include.lowest = TRUE, labels = FALSE)  
  return(list(y = y, X = X, k=(length(cutoff)-1), betas= betas, cutoff=cutoff ))
}

nobs = 300 
X = cbind(rep(1,nobs),runif(nobs, min=0, max=5),runif(nobs,min=0, max=5))
k = 5
betas = c(0.5, 1, -0.5)       
cutoff = c(-100, 0, 1.0, 1.8, 3.2, 100)
simout = simordprobit(X, betas, cutoff) 
  
Data=list(X=simout$X, y=simout$y, k=k)

## set Mcmc for ordered probit model
   
Mcmc = list(R=R)   
out = rordprobitGibbs(Data=Data, Mcmc=Mcmc)

cat(" ", fill=TRUE)
cat("acceptance rate= ", accept=out$accept, fill=TRUE)
 
## outputs of betadraw and cut-off draws
  
cat(" Summary of betadraws", fill=TRUE)
summary(out$betadraw, tvalues=betas)
cat(" Summary of cut-off draws", fill=TRUE) 
summary(out$cutdraw, tvalues=cutoff[2:k])

## plotting examples
if(0){plot(out$cutdraw)}
}

\keyword{models}
