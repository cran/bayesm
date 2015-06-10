\name{rnegbinRw}
\alias{rnegbinRw}
\concept{MCMC}
\concept{NBD regression}
\concept{Negative Binomial regression}
\concept{Poisson regression}
\concept{Metropolis algorithm}
\concept{bayes}
\title{ MCMC Algorithm for Negative Binomial Regression }
\description{
 \code{rnegbinRw} implements a Random Walk Metropolis Algorithm for the Negative 
 Binomial (NBD) regression model. beta | alpha and alpha | beta are drawn with two
 different random walks.
}
\usage{
rnegbinRw(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y,X) }
  \item{Prior}{ list(betabar,A,a,b) }
  \item{Mcmc}{ list(R,keep,s_beta,s_alpha,beta0 }
}
\details{
  Model:   \eqn{y} \eqn{\sim}{~} \eqn{NBD(mean=\lambda, over-dispersion=alpha)}.  \cr
           \eqn{\lambda=exp(x'\beta)}

  Prior:   \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar,A^{-1})} \cr
           \eqn{alpha} \eqn{\sim}{~} \eqn{Gamma(a,b)}. \cr
            note: prior mean of \eqn{alpha = a/b}, \eqn{variance = a/(b^2)}

  list arguments contain:
  \itemize{
    \item{\code{y}}{ nobs vector of counts (0,1,2,\ldots)}
    \item{\code{X}}{nobs x nvar matrix}
    \item{\code{betabar}}{ nvar x 1 prior mean (def: 0)}
    \item{\code{A}}{ nvar x nvar pds prior prec matrix (def: .01I)}
    \item{\code{a}}{ Gamma prior parm (def: .5)}
    \item{\code{b}}{ Gamma prior parm (def: .1)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
    \item{\code{s_beta}}{ scaling for beta| alpha RW inc cov matrix (def: 2.93/sqrt(nvar))}
    \item{\code{s_alpha}}{ scaling for alpha | beta RW inc cov matrix (def: 2.93)}
  }
}
\value{
  a list containing: 
  \item{betadraw}{R/keep x nvar array of beta draws}
  \item{alphadraw}{R/keep vector of alpha draws}
  \item{llike}{R/keep vector of log-likelihood values evaluated at each draw}
  \item{acceptrbeta}{acceptance rate of the beta draws}
  \item{acceptralpha}{acceptance rate of the alpha draws}
}
\note{
  The NBD regression encompasses Poisson regression in the sense that as alpha goes to
  infinity the NBD distribution tends toward the Poisson.\cr
  For "small" values of alpha, the dependent variable can be extremely variable so that 
  a large number of observations may be required to obtain precise inferences.
}

\seealso{ \code{\link{rhierNegbinRw}} }

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby, McCulloch. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Sridhar Narayanam & Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0)  {R=1000} else {R=10}

set.seed(66)
simnegbin = 
function(X, beta, alpha) {
#   Simulate from the Negative Binomial Regression
lambda = exp(X \%*\% beta)
y=NULL
for (j in 1:length(lambda))
    y = c(y,rnbinom(1,mu = lambda[j],size = alpha))
return(y)
}

nobs = 500
nvar=2            # Number of X variables
alpha = 5
Vbeta = diag(nvar)*0.01

# Construct the regdata (containing X)
simnegbindata = NULL
beta = c(0.6,0.2)
X = cbind(rep(1,nobs),rnorm(nobs,mean=2,sd=0.5))
simnegbindata = list(y=simnegbin(X,beta,alpha), X=X, beta=beta)

Data1 = simnegbindata
Mcmc1 = list(R=R)

out = rnegbinRw(Data=Data1,Mcmc=Mcmc1)

cat("Summary of alpha/beta draw",fill=TRUE)
summary(out$alphadraw,tvalues=alpha)
summary(out$betadraw,tvalues=beta)

if(0){
## plotting examples
plot(out$betadraw)
}

}

\keyword{ models }

