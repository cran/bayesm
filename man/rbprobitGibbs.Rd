\name{rbprobitGibbs}
\alias{rbprobitGibbs}
\concept{bayes}
\concept{MCMC}
\concept{probit}
\concept{Gibbs Sampling}

\title{ Gibbs Sampler (Albert and Chib) for Binary Probit }
\description{
  \code{rbprobitGibbs} implements the Albert and Chib Gibbs Sampler for the binary probit model.

}
\usage{
rbprobitGibbs(Data, Prior, Mcmc)
}

\arguments{
  \item{Data}{ list(X,y)}
  \item{Prior}{ list(betabar,A)}
  \item{Mcmc}{ list(R,keep)  }
}

\details{
  Model: \eqn{z = X\beta + e}.  \eqn{e} \eqn{\sim}{~} \eqn{N(0,I)}. y=1, if z> 0. 

  Prior:  \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})}. 

  List arguments contain  
  \describe{
    \item{\code{X}}{Design Matrix}
    \item{\code{y}}{n x 1 vector of observations, (0 or 1)} 
    \item{\code{betabar}}{k x 1 prior mean (def: 0)}
    \item{\code{A}}{k x k prior precision matrix (def: .01I)} 
    \item{\code{R}}{ number of MCMC draws }
    \item{\code{keep}}{ thinning parameter - keep every keepth draw (def: 1)}
  }
}

\value{
  \item{betadraw }{R/keep x k array of betadraws}
}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 3. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\seealso{ \code{\link{rmnpGibbs}} }
\examples{
##
## rbprobitGibbs example
##
set.seed(66)
simbprobit=
function(X,beta) {
##  function to simulate from binary probit including x variable
y=ifelse((X\%*\%beta+rnorm(nrow(X)))<0,0,1)
list(X=X,y=y,beta=beta)
}

nobs=200
X=cbind(rep(1,nobs),runif(nobs),runif(nobs))
beta=c(0,1,-1)
nvar=ncol(X)
simout=simbprobit(X,beta)

Data=list(X=simout$X,y=simout$y)
Mcmc=list(R=2000,keep=1)

out=rbprobitGibbs(Data=Data,Mcmc=Mcmc)

cat(" Betadraws ",fill=TRUE)
mat=apply(out$betadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(beta,mat); rownames(mat)[1]="beta"; print(mat)
}
\keyword{ models }
