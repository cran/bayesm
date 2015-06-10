\name{rsurGibbs}
\alias{rsurGibbs}
\concept{bayes}
\concept{Gibbs Sampler}
\concept{regression}
\concept{SUR model}
\concept{Seemingly Unrelated Regression}
\concept{MCMC}

\title{ Gibbs Sampler for Seemingly Unrelated Regressions (SUR) }
\description{
 \code{rsurGibbs} implements a Gibbs Sampler to draw from the posterior of the Seemingly Unrelated
  Regression (SUR) Model of Zellner
}
\usage{
rsurGibbs(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(regdata)}
  \item{Prior}{ list(betabar,A, nu, V) }
  \item{Mcmc}{ list(R,keep)}
}
\details{
  Model: \eqn{y_i = X_i\beta_i + e_i}.  i=1,\ldots,m. m regressions. \cr
  (e(1,k), \ldots, e(m,k)) \eqn{\sim}{~} \eqn{N(0,\Sigma)}. k=1, \ldots, nobs. 

  We can also write as the stacked model: \cr
  \eqn{y = X\beta + e} where y is a nobs*m long vector and k=length(beta)=sum(length(betai)).

  Note: we must have the same number of observations in each equation but we can have different numbers
  of X variables 

  Priors: \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})}.  \eqn{\Sigma} \eqn{\sim}{~} \eqn{IW(nu,V)}.

  List arguments contain  
  \itemize{
    \item{\code{regdata}}{list of lists, regdata[[i]]=list(y=yi,X=Xi)}
    \item{\code{betabar}}{k x 1 prior mean (def: 0)}
    \item{\code{A}}{k x k prior precision matrix (def: .01I)} 
    \item{\code{nu}}{ d.f. parm for Inverted Wishart prior (def: m+3)}
    \item{\code{V}}{ scale parm for Inverted Wishart prior (def: nu*I)}
    \item{\code{R}}{ number of MCMC draws }
    \item{\code{keep}}{ thinning parameter - keep every keepth draw }
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
  }
}
\value{
  list of MCMC draws
  \item{betadraw }{ R x k array of betadraws }
  \item{Sigmadraw }{ R x (m*m) array of Sigma draws}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 3. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\seealso{ \code{\link{rmultireg}} }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=1000} else {R=10}
##
## simulate data from SUR
set.seed(66)
beta1=c(1,2)
beta2=c(1,-1,-2)
nobs=100
nreg=2
iota=c(rep(1,nobs))
X1=cbind(iota,runif(nobs))
X2=cbind(iota,runif(nobs),runif(nobs))
Sigma=matrix(c(.5,.2,.2,.5),ncol=2)
U=chol(Sigma)
E=matrix(rnorm(2*nobs),ncol=2)\%*\%U
y1=X1\%*\%beta1+E[,1]
y2=X2\%*\%beta2+E[,2]
##
## run Gibbs Sampler
regdata=NULL
regdata[[1]]=list(y=y1,X=X1)
regdata[[2]]=list(y=y2,X=X2)

Mcmc1=list(R=R)

out=rsurGibbs(Data=list(regdata=regdata),Mcmc=Mcmc1)

cat("Summary of beta draws",fill=TRUE)
summary(out$betadraw,tvalues=c(beta1,beta2))
cat("Summary of Sigmadraws",fill=TRUE)
summary(out$Sigmadraw,tvalues=as.vector(Sigma[upper.tri(Sigma,diag=TRUE)]))

if(0){
plot(out$betadraw,tvalues=c(beta1,beta2))
}

}

\keyword{ regression}
