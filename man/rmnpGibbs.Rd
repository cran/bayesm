\name{rmnpGibbs}
\alias{rmnpGibbs}
\concept{bayes}
\concept{multinomial probit}
\concept{MCMC}
\concept{Gibbs Sampling}

\title{ Gibbs Sampler for Multinomial Probit }
\description{
  \code{rmnpGibbs} implements the McCulloch/Rossi Gibbs Sampler for the multinomial probit model.
}

\usage{
rmnpGibbs(Data, Prior, Mcmc)
}

\arguments{
  \item{Data}{ list(p, y, X)}
  \item{Prior}{ list(betabar,A,nu,V) (optional)}
  \item{Mcmc}{ list(beta0,sigma0,R,keep,nprint) (R required) }
}

\details{
  model:  \cr
    \eqn{w_i = X_i\beta + e}.    \eqn{e} \eqn{\sim}{~} \eqn{N(0,\Sigma)}.     note: \eqn{w_i, e} are (p-1) x 1.\cr
    \eqn{y_i = j},  if \eqn{w_{ij} > max(0,w_{i,-j})}  j=1,\ldots,p-1.  \eqn{w_{i,-j}} means elements of \eqn{w_i}
     other than the jth. \cr
    \eqn{y_i = p},  if all \eqn{w_i < 0}.\cr
  
  priors:\cr
    \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})} \cr
    \eqn{\Sigma} \eqn{\sim}{~} IW(nu,V)\cr

  to make up X matrix use \code{\link{createX}} with \code{DIFF=TRUE}.

  List arguments contain  
  \itemize{
    \item{\code{p}}{number of choices or possible multinomial outcomes}
    \item{\code{y}}{n x 1 vector of multinomial outcomes}
    \item{\code{X}}{n*(p-1) x k Design Matrix}
    \item{\code{betabar}}{k x 1 prior mean (def: 0)}
    \item{\code{A}}{k x k prior precision matrix (def: .01I)} 
    \item{\code{nu}}{ d.f. parm for IWishart prior (def: (p-1) + 3)}
    \item{\code{V}}{ pds location parm for IWishart prior (def: nu*I)}
    \item{\code{beta0}}{ initial value for beta}
    \item{\code{sigma0}}{ initial value for sigma }
    \item{\code{R}}{ number of MCMC draws }
    \item{\code{keep}}{ thinning parameter - keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
  }
}

\value{
  a list containing: 
  \item{betadraw }{R/keep x k array of betadraws}
  \item{sigmadraw}{R/keep x (p-1)*(p-1) array of sigma draws -- each row is in vector form}
}
\note{
  \eqn{\beta} is not identified.  \eqn{\beta}/sqrt(\eqn{\sigma_{11}}) and \eqn{\Sigma}/\eqn{\sigma_{11}} are.  See Allenby et al or
  example below for details.
}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 4. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\seealso{ \code{\link{rmvpGibbs}} }
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(66)
p=3
n=500
beta=c(-1,1,1,2)
Sigma=matrix(c(1,.5,.5,1),ncol=2)
k=length(beta)
X1=matrix(runif(n*p,min=0,max=2),ncol=p); X2=matrix(runif(n*p,min=0,max=2),ncol=p)
X=createX(p,na=2,nd=NULL,Xa=cbind(X1,X2),Xd=NULL,DIFF=TRUE,base=p)

simmnp= function(X,p,n,beta,sigma) {
  indmax=function(x) {which(max(x)==x)}
  Xbeta=X\%*\%beta
  w=as.vector(crossprod(chol(sigma),matrix(rnorm((p-1)*n),ncol=n)))+ Xbeta
  w=matrix(w,ncol=(p-1),byrow=TRUE)
  maxw=apply(w,1,max)
  y=apply(w,1,indmax)
  y=ifelse(maxw < 0,p,y)
  return(list(y=y,X=X,beta=beta,sigma=sigma))
}

simout=simmnp(X,p,500,beta,Sigma)

Data1=list(p=p,y=simout$y,X=simout$X)
Mcmc1=list(R=R,keep=1)

out=rmnpGibbs(Data=Data1,Mcmc=Mcmc1)

cat(" Summary of Betadraws ",fill=TRUE)
betatilde=out$betadraw/sqrt(out$sigmadraw[,1])
attributes(betatilde)$class="bayesm.mat"
summary(betatilde,tvalues=beta)

cat(" Summary of Sigmadraws ",fill=TRUE)
sigmadraw=out$sigmadraw/out$sigmadraw[,1]
attributes(sigmadraw)$class="bayesm.var"
summary(sigmadraw,tvalues=as.vector(Sigma[upper.tri(Sigma,diag=TRUE)]))


if(0){
## plotting examples
plot(betatilde,tvalues=beta)
}
}
\keyword{ models }
