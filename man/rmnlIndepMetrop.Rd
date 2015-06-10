\name{rmnlIndepMetrop}
\alias{rmnlIndepMetrop}
\concept{MCMC}
\concept{multinomial logit}
\concept{Metropolis algorithm}
\concept{bayes}
\title{ MCMC Algorithm for Multinomial Logit Model }
\description{
  \code{rmnIndepMetrop} implements Independence Metropolis for the MNL.
}
\usage{
rmnlIndepMetrop(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(p,y,X)}
  \item{Prior}{ list(A,betabar)  optional}
  \item{Mcmc}{ list(R,keep,nprint,nu) }
}
\details{
  Model:   y \eqn{\sim}{~} MNL(X,\eqn{\beta}). \eqn{\Pr(y=j) = exp(x_j'\beta)/\sum_k{e^{x_k'\beta}}}. \cr

  Prior:   \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})} \cr

  list arguments contain:
  \itemize{
    \item{\code{p}}{ number of alternatives}
    \item{\code{y}}{ nobs vector of multinomial outcomes (1,\ldots, p)}
    \item{\code{X}}{ nobs*p x nvar matrix}
    \item{\code{A}}{ nvar x nvar pds prior prec matrix (def: .01I)}
    \item{\code{betabar}}{ nvar x 1 prior mean (def: 0)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
    \item{\code{nu}}{ degrees of freedom parameter for independence t density (def: 6) }
  }
}
\value{
  a list containing: 
  \item{betadraw}{R/keep x nvar array of beta draws}
  \item{loglike}{R/keep vector of loglike values for each draw}
  \item{acceptr}{acceptance rate of Metropolis draws}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 5. \cr
  \url{http://www.perossi.org/home/bsm-1}
}
\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\seealso{ \code{\link{rhierMnlRwMixture}} }
\examples{
##

if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(66)
n=200; p=3; beta=c(1,-1,1.5,.5)

simmnl= function(p,n,beta) {
  #   note: create X array with 2 alt.spec vars
  k=length(beta)
  X1=matrix(runif(n*p,min=-1,max=1),ncol=p)
  X2=matrix(runif(n*p,min=-1,max=1),ncol=p)
  X=createX(p,na=2,nd=NULL,Xd=NULL,Xa=cbind(X1,X2),base=1)
  Xbeta=X\%*\%beta # now do probs
  p=nrow(Xbeta)/n
  Xbeta=matrix(Xbeta,byrow=TRUE,ncol=p)
  Prob=exp(Xbeta)
  iota=c(rep(1,p))
  denom=Prob\%*\%iota
  Prob=Prob/as.vector(denom)
  # draw y
  y=vector("double",n)
  ind=1:p
  for (i in 1:n) 
        { yvec=rmultinom(1,1,Prob[i,]); y[i]=ind\%*\%yvec }
   return(list(y=y,X=X,beta=beta,prob=Prob))
}

simout=simmnl(p,n,beta)

Data1=list(y=simout$y,X=simout$X,p=p); Mcmc1=list(R=R,keep=1)
out=rmnlIndepMetrop(Data=Data1,Mcmc=Mcmc1)

cat("Summary of beta draws",fill=TRUE)
summary(out$betadraw,tvalues=beta)

if(0){
## plotting examples
plot(out$betadraw)
}

}
\keyword{ models }
