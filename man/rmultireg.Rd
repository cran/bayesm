\name{rmultireg}
\alias{rmultireg}
\concept{bayes}
\concept{multivariate regression}
\concept{simulation}


\title{ Draw from the Posterior of a Multivariate Regression }
\description{
  \code{ rmultireg} draws from the posterior of a Multivariate Regression model with a natural conjugate
  prior.
}
\usage{
rmultireg(Y, X, Bbar, A, nu, V)
}

\arguments{
  \item{Y}{ n x m matrix of observations on m dep vars }
  \item{X}{ n x k matrix of observations on indep vars (supply intercept) }
  \item{Bbar}{ k x m matrix of prior mean of regression coefficients }
  \item{A}{ k x k Prior precision matrix }
  \item{nu}{ d.f. parameter for Sigma }
  \item{V}{ m x m pdf location parameter for prior on Sigma }
}
\details{
  Model: \eqn{Y=XB+U}.  \eqn{cov(u_i) = \Sigma}.  \eqn{B} is k x m matrix of coefficients. \eqn{\Sigma} is m x m covariance.

  Priors:  \eqn{\beta} given \eqn{\Sigma}  \eqn{\sim}{~} \eqn{N(betabar,\Sigma (x) A^{-1})}. 
  \eqn{betabar=vec(Bbar)};  \eqn{\beta = vec(B)} \cr
          \eqn{\Sigma} \eqn{\sim}{~} IW(nu,V). 
}
\value{
  A list of the components of a draw from the posterior
  \item{B }{ draw of regression coefficient matrix }
  \item{Sigma }{ draw of Sigma }
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 2. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(66)
n=200
m=2
X=cbind(rep(1,n),runif(n))
k=ncol(X)
B=matrix(c(1,2,-1,3),ncol=m)
Sigma=matrix(c(1,.5,.5,1),ncol=m); RSigma=chol(Sigma)
Y=X\%*\%B+matrix(rnorm(m*n),ncol=m)\%*\%RSigma

betabar=rep(0,k*m);Bbar=matrix(betabar,ncol=m)
A=diag(rep(.01,k))
nu=3; V=nu*diag(m)

betadraw=matrix(double(R*k*m),ncol=k*m)
Sigmadraw=matrix(double(R*m*m),ncol=m*m)
for (rep in 1:R) 
   {out=rmultireg(Y,X,Bbar,A,nu,V);betadraw[rep,]=out$B
    Sigmadraw[rep,]=out$Sigma}

cat(" Betadraws ",fill=TRUE)
mat=apply(betadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(as.vector(B),mat); rownames(mat)[1]="beta"
print(mat)
cat(" Sigma draws",fill=TRUE)
mat=apply(Sigmadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(as.vector(Sigma),mat); rownames(mat)[1]="Sigma"
print(mat)
}
\keyword{ regression }
