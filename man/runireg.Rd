\name{runireg}
\alias{runireg}
\concept{bayes}
\concept{regression}
\concept{simulation}

\title{ Draw from Posterior for Univariate Regression }
\description{
 \code{runireg} draws from posterior for univariate regression with a natural conjugate prior.
}
\usage{
runireg(y, X, betabar, A, nu, ssq)
}
\arguments{
  \item{y}{ n x 1 dep var  }
  \item{X}{ n x k Design matrix }
  \item{betabar}{ prior mean }
  \item{A}{ k x k pds prior precision matrix }
  \item{nu}{ d.f. parameter for sigma-sq prior }
  \item{ssq}{ scale parameter for sigma-sq prior }
}
\details{
  Model: \eqn{y = Xbeta + e}.  \eqn{e} \eqn{\sim}{~} \eqn{N(0,sigmasq)}. \cr

  Priors: beta given sigmasq \eqn{\sim}{~} \eqn{N(betabar,sigmasq*A^{-1})}.\cr
  \eqn{sigmasq} \eqn{\sim}{~} eqn{(nu*ssq)}/\eqn{\chi^2_{nu}}.
}
\value{
  a list with one draw from posterior
  \item{beta }{beta draw}
  \item{sigmasq }{sigmasq draw}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 2. \cr 
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{runiregGibbs}} }
\examples{
set.seed(66)
n=100
X=cbind(rep(1,n),runif(n)); beta=c(1,2); sigsq=.25
y=X\%*\%beta+rnorm(n,sd=sqrt(sigsq))

A=diag(c(.05,.05)); betabar=c(0,0)
nu=3; ssq=1.0

R=1000
betadraw=matrix(double(R*2),ncol=2)
sigsqdraw=double(R)
for (rep in 1:R) 
   {out=runireg(y,X,betabar,A,nu,ssq);betadraw[rep,]=out$beta
    sigsqdraw[rep]=out$sigmasq}

cat(" Betadraws ",fill=TRUE)
mat=apply(betadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(beta,mat); rownames(mat)[1]="beta"; print(mat)
cat(" Sigma-sq draws",fill=TRUE)
cat(" sigma-sq= ",sigsq,fill=TRUE)
print(quantile(sigsqdraw,probs=c(.01,.05,.5,.95,.99)))
}
\keyword{ regression }
