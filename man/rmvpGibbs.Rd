\name{rmvpGibbs}
\alias{rmvpGibbs}

\concept{bayes}
\concept{multivariate probit}
\concept{MCMC}
\concept{Gibbs Sampling}

\title{ Gibbs Sampler for Multivariate Probit }
\description{
  \code{rmvpGibbs} implements the Edwards/Allenby Gibbs Sampler for the multivariate probit model.
}
  
}
\usage{
rmvpGibbs(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(p,y,X)}
  \item{Prior}{ list(betabar,A,nu,V) (optional)}
  \item{Mcmc}{ list(beta0,sigma0,R,keep) (R required) }
}

\details{
  model:  \cr
    \eqn{w_i = X_i beta + e}.    \eqn{e} \eqn{\sim}{~} N(0,Sigma).     note: \eqn{w_i} is p x 1.\cr
    \eqn{y_{ij} = 1},  if \eqn{w_{ij} > 0}, else \eqn{y_i=0}.  j=1,\ldots,p.   \cr
  
  priors:\cr
    \eqn{beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})}\cr
    \eqn{Sigma} \eqn{\sim}{~} IW(nu,V)\cr

  to make up X matrix use \code{createX}

  List arguments contain  
  \itemize{
    \item{\code{p}}{dimension of multivariate probit}
    \item{\code{X}}{n*p x k Design Matrix}
    \item{\code{y}}{n*p x 1 vector of 0,1 outcomes}
    \item{\code{betabar}}{k x 1 prior mean (def: 0)}
    \item{\code{A}}{k x k prior precision matrix (def: .01I)} 
    \item{\code{nu}}{ d.f. parm for IWishart prior (def: (p-1) + 3)}
    \item{\code{V}}{ pds location parm for IWishart prior (def: nu*I)}
    \item{\code{beta0}}{ initial value for beta}
    \item{\code{sigma0}}{ initial value for sigma }
    \item{\code{R}}{ number of MCMC draws }
    \item{\code{keep}}{ thinning parameter - keep every keepth draw (def: 1)}
  }
}

\value{
  a list containing: 
  \item{betadraw }{R/keep x k array of betadraws}
  \item{sigmadraw}{R/keep x p*p array of sigma draws -- each row is in vector form}
}

\note{
  beta and Sigma are not identifed.  Correlation matrix and the betas divided by the
  appropriate standard deviation are.  See Allenby et al for details or example below.
} 

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 4. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\seealso{ \code{\link{rmnpGibbs}} }
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(66)
p=3
n=500
beta=c(-2,0,2)
Sigma=matrix(c(1,.5,.5,.5,1,.5,.5,.5,1),ncol=3)
k=length(beta)
I2=diag(rep(1,p)); xadd=rbind(I2)
for(i in 2:n) { xadd=rbind(xadd,I2)}; X=xadd
simout=simmvp(X,p,500,beta,Sigma)

Data=list(p=p,y=simout$y,X=simout$X)
Mcmc=list(R=R,keep=1)
out=rmvpGibbs(Data=Data,Mcmc=Mcmc)

ind=seq(from=0,by=p,length=k)
inda=1:3
ind=ind+inda
cat(" Betadraws ",fill=TRUE)
mat=apply(out$betadraw/sqrt(out$sigmadraw[,ind]),2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(beta,mat); rownames(mat)[1]="beta"; print(mat)
rdraw=matrix(double((R)*p*p),ncol=p*p)
rdraw=t(apply(out$sigmadraw,1,nmat))
cat(" Draws of Correlation Matrix ",fill=TRUE)
mat=apply(rdraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(as.vector(Sigma),mat); rownames(mat)[1]="Sigma"; print(mat)

}
\keyword{ models }
\keyword{ multivariate }
