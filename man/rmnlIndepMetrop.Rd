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
  \item{Mcmc}{ list(R,keep,nu) }
}
\details{
  Model:   y \eqn{\sim}{~} MNL(X,beta). \eqn{Pr(y=j) = exp(x_j'beta)/\sum_k{e^{x_k'beta}}}. \cr

  Prior:   \eqn{beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})} \cr

  list arguments contain:
  \itemize{
    \item{\code{p}}{number of alternatives}
    \item{\code{y}}{ nobs vector of multinomial outcomes (1,\ldots, p)}
    \item{\code{X}}{nobs*p x nvar matrix}
    \item{\code{A}}{ nvar x nvar pds prior prec matrix (def: .01I)}
    \item{\code{betabar}}{ nvar x 1 prior mean (def: 0)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
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
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}
\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}
\seealso{ \code{\link{rhierMnlRwMixture}} }
\examples{
##

if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(66)
n=200; p=3; beta=c(1,-1,1.5,.5)
simout=simmnl(p,n,beta)
A=diag(c(rep(.01,length(beta)))); betabar=rep(0,length(beta))

Data=list(y=simout$y,X=simout$X,p=p); Mcmc=list(R=R,keep=1) ; Prior=list(A=A,betabar=betabar)
out=rmnlIndepMetrop(Data=Data,Prior=Prior,Mcmc=Mcmc)
cat(" Betadraws ",fill=TRUE)
mat=apply(out$betadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(beta,mat); rownames(mat)[1]="beta"; print(mat)

}
\keyword{ models }
