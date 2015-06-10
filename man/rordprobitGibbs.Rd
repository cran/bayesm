\name{rordprobitGibbs}
\alias{rordprobitGibbs}
\concept{bayes}
\concept{MCMC}
\concept{probit}
\concept{Gibbs Sampling}

\title{ Gibbs Sampler for Ordered Probit }
\description{
  \code{rordprobitGibbs} implements a Gibbs Sampler for the ordered probit model.

}
\usage{
rordprobitGibbs(Data, Prior, Mcmc)
}

\arguments{
  \item{Data}{ list(X, y, k)}
  \item{Prior}{ list(betabar, A, dstarbar, Ad)}
  \item{Mcmc}{ list(R, keep, nprint s, change, draw)  }
}

\details{
  Model: \eqn{z = X\beta + e}.  \eqn{e} \eqn{\sim}{~} \eqn{N(0,I)}. 
          y=1,..,k. cutoff=c( c [1] ,..c [k+1] ).   \cr 
          y=k, if c [k] <= z < c [k+1] .     

  Prior:  \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})}. 
          \eqn{dstar} \eqn{\sim}{~} \eqn{N(dstarbar,Ad^{-1})}.

  List arguments contain  
  \describe{
    \item{\code{X}}{n x nvar Design Matrix}
    \item{\code{y}}{n x 1 vector of observations, (1,...,k)} 
    \item{\code{k}}{the largest possible value of y} 
    \item{\code{betabar}}{nvar x 1 prior mean (def: 0)}
    \item{\code{A}}{nvar x nvar prior precision matrix (def: .01I)} 
    \item{\code{dstarbar}}{ndstar x 1 prior mean, ndstar=k-2 (def: 0)}
    \item{\code{Ad}}{ndstar x ndstar prior precision matrix (def:I)} 
    \item{\code{s}}{ scaling parm for RW Metropolis (def: 2.93/sqrt(nvar))}
    \item{\code{R}}{ number of MCMC draws }
    \item{\code{keep}}{ thinning parameter - keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
  }
}

\value{
  \item{betadraw }{R/keep x k matrix of betadraws}
  \item{cutdraw }{R/keep x (k-1) matrix of cutdraws}
  \item{dstardraw }{R/keep x (k-2) matrix of dstardraws}
  \item{accept }{a value of acceptance rate in RW Metropolis}
}
\note{ 
   set c[1]=-100. c[k+1]=100. c[2] is set to 0 for identification.   \cr

   The relationship between cut-offs and dstar is    \cr
   c[3] = exp(dstar[1]), c[4]=c[3]+exp(dstar[2]),..., c[k] = c[k-1] + exp(datsr[k-2])              
  
   Be careful in assessing prior parameter, Ad.  .1 is too small for many applications. 
} 

\references{  \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch\cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\seealso{ \code{\link{rbprobitGibbs}} }
\examples{
##
## rordprobitGibbs example
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

## simulate data for ordered probit model

   simordprobit=function(X, betas, cutoff){
    z = X\%*\%betas + rnorm(nobs)   
    y = cut(z, br = cutoff, right=TRUE, include.lowest = TRUE, labels = FALSE)  
    return(list(y = y, X = X, k=(length(cutoff)-1), betas= betas, cutoff=cutoff ))
   }

   set.seed(66)  
   nobs=300 
   X=cbind(rep(1,nobs),runif(nobs, min=0, max=5),runif(nobs,min=0, max=5))
   k=5
   betas=c(0.5, 1, -0.5)       
   cutoff=c(-100, 0, 1.0, 1.8, 3.2,  100)
   simout=simordprobit(X, betas, cutoff)   
   Data=list(X=simout$X,y=simout$y, k=k)

## set Mcmc for ordered probit model
   
   Mcmc=list(R=R)   
   out=rordprobitGibbs(Data=Data,Mcmc=Mcmc)
  
   cat(" ", fill=TRUE)
   cat("acceptance rate= ",accept=out$accept,fill=TRUE)
 
## outputs of betadraw and cut-off draws
  
   cat(" Summary of betadraws",fill=TRUE)
   summary(out$betadraw,tvalues=betas)
   cat(" Summary of cut-off draws",fill=TRUE) 
   summary(out$cutdraw,tvalues=cutoff[2:k])

if(0){
## plotting examples
plot(out$cutdraw)
}

}
\keyword{ models }
