\name{rhierLinearModel}
\alias{rhierLinearModel}
\concept{bayes}
\concept{MCMC}
\concept{Gibbs Sampling}
\concept{hierarchical models}
\concept{linear model}

\title{ Gibbs Sampler for Hierarchical Linear Model }
\description{
  \code{rhierLinearModel} implements a Gibbs Sampler for hierarchical linear models with a normal prior.
}
\usage{
rhierLinearModel(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(regdata,Z) (Z optional). }
  \item{Prior}{ list(Deltabar,A,nu.e,ssq,nu,V)  (optional).}
  \item{Mcmc}{ list(R,keep,nprint) (R required).}
}
\details{
  Model: length(regdata) regression equations. \cr
        \eqn{y_i = X_i\beta_i + e_i}. \eqn{e_i} \eqn{\sim}{~} \eqn{N(0,\tau_i)}.  nvar X vars in each equation. 

 Priors:\cr
        \eqn{\tau_i} \eqn{\sim}{~} nu.e*\eqn{ssq_i/\chi^2_{nu.e}}.  \eqn{\tau_i} is the variance of \eqn{e_i}.\cr
        \eqn{\beta_i} \eqn{\sim}{~} N(Z\eqn{\Delta}[i,],\eqn{V_{\beta}}). \cr
               Note:  Z\eqn{\Delta} is the matrix Z * \eqn{\Delta}; [i,] refers to ith row of this product.\cr

          \eqn{vec(\Delta)} given \eqn{V_{\beta}} \eqn{\sim}{~} \eqn{N(vec(Deltabar),V_{\beta} (x) A^{-1})}.\cr
          \eqn{V_{\beta}} \eqn{\sim}{~} \eqn{IW(nu,V)}. \cr
              \eqn{Delta, Deltabar} are nz x nvar.  \eqn{A} is nz x nz.  \eqn{V_{\beta}} is nvar x nvar.
        
          Note: if you don't have any Z vars, omit Z in the \code{Data} argument and a vector of ones will be inserted for you.  In this case (of no Z vars), the matrix \eqn{\Delta} will be 1 x nvar and should be interpreted as the mean of all unit \eqn{\beta} s.

  List arguments contain:
  \itemize{
    \item{\code{regdata}}{ list of lists with X,y matrices for each of length(regdata) regressions}
    \item{\code{regdata[[i]]$X}}{ X matrix for equation i }
    \item{\code{regdata[[i]]$y}}{ y vector for equation i }
    \item{\code{Deltabar}}{ nz x nvar matrix of prior means (def: 0)}
    \item{\code{A}}{ nz x nz matrix for prior precision (def: .01I)}
    \item{\code{nu.e}}{ d.f. parm for regression error variance prior (def: 3)}
    \item{\code{ssq}}{ scale parm for regression error var prior (def: var(\eqn{y_i}))}
    \item{\code{nu}}{ d.f. parm for Vbeta prior (def: nvar+3)}
    \item{\code{V}}{ Scale location matrix for Vbeta prior (def: nu*I)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
   }
}
\value{
  a list containing
  \item{betadraw}{nreg x nvar x R/keep array of individual regression coef draws}
  \item{taudraw}{R/keep x nreg array of error variance draws}
  \item{Deltadraw}{R/keep x nz x nvar array of Deltadraws}
  \item{Vbetadraw}{R/keep x nvar*nvar array of Vbeta draws}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 3. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.comu}.
}
\seealso{ \code{\link{rhierLinearMixture}} }
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

nreg=100; nobs=100; nvar=3
Vbeta=matrix(c(1,.5,0,.5,2,.7,0,.7,1),ncol=3)
Z=cbind(c(rep(1,nreg)),3*runif(nreg)); Z[,2]=Z[,2]-mean(Z[,2])
nz=ncol(Z)
Delta=matrix(c(1,-1,2,0,1,0),ncol=2)
Delta=t(Delta) # first row of Delta is means of betas
Beta=matrix(rnorm(nreg*nvar),nrow=nreg)\%*\%chol(Vbeta)+Z\%*\%Delta
tau=.1
iota=c(rep(1,nobs))
regdata=NULL
for (reg in 1:nreg) { X=cbind(iota,matrix(runif(nobs*(nvar-1)),ncol=(nvar-1)))
	y=X\%*\%Beta[reg,]+sqrt(tau)*rnorm(nobs); regdata[[reg]]=list(y=y,X=X) }

Data1=list(regdata=regdata,Z=Z)
Mcmc1=list(R=R,keep=1)
out=rhierLinearModel(Data=Data1,Mcmc=Mcmc1)

cat("Summary of Delta draws",fill=TRUE)
summary(out$Deltadraw,tvalues=as.vector(Delta))
cat("Summary of Vbeta draws",fill=TRUE)
summary(out$Vbetadraw,tvalues=as.vector(Vbeta[upper.tri(Vbeta,diag=TRUE)]))

if(0){
## plotting examples
plot(out$betadraw)
plot(out$Deltadraw)
}

}
\keyword{ regression }
