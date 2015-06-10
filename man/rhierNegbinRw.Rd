\name{rhierNegbinRw}
\alias{rhierNegbinRw}
\concept{MCMC}
\concept{hierarchical NBD regression}
\concept{Negative Binomial regression}
\concept{Poisson regression}
\concept{Metropolis algorithm}
\concept{bayes}
\title{ MCMC Algorithm for Negative Binomial Regression }
\description{
 \code{rhierNegbinRw} implements an MCMC strategy for the hierarchical Negative 
 Binomial (NBD) regression model. Metropolis steps for each unit level set of 
 regression parameters are automatically tuned by optimization. Over-dispersion
 parameter (alpha) is common across units.
}
\usage{
rhierNegbinRw(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(regdata,Z) }
  \item{Prior}{ list(Deltabar,Adelta,nu,V,a,b) }
  \item{Mcmc}{ list(R,keep,nprint,s_beta,s_alpha,c,Vbeta0,Delta0) }
}
\details{
  Model:   \eqn{y_i} \eqn{\sim}{~} NBD(mean=\eqn{\lambda}, over-dispersion=alpha).  \cr
           
           \eqn{\lambda=exp(X_i\beta_i)}

  Prior:   \eqn{\beta_i} \eqn{\sim}{~} \eqn{N(\Delta'z_i,Vbeta)}.

           \eqn{vec(\Delta|Vbeta)} \eqn{\sim}{~} \eqn{N(vec(Deltabar),Vbeta (x) Adelta)}. \cr
           \eqn{Vbeta} \eqn{\sim}{~} \eqn{IW(nu,V)}. \cr
           \eqn{alpha} \eqn{\sim}{~} \eqn{Gamma(a,b)}. \cr
            note: prior mean of \eqn{alpha = a/b}, \eqn{variance = a/(b^2)}

  list arguments contain:
  \itemize{
    \item{\code{regdata}}{ list of lists with data on each of nreg units}
    \item{\code{regdata[[i]]$X}}{ nobs\_i x nvar matrix of X variables}
    \item{\code{regdata[[i]]$y}}{ nobs\_i x 1 vector of count responses}
    \item{\code{Z}}{nreg x nz mat of unit chars (def: vector of ones)}
    \item{\code{Deltabar}}{ nz x nvar prior mean matrix (def: 0)}
    \item{\code{Adelta}}{ nz x nz pds prior prec matrix (def: .01I)}
    \item{\code{nu}}{ d.f. parm for IWishart (def: nvar+3)}
    \item{\code{V}}{location matrix of IWishart prior (def: nuI)}
    \item{\code{a}}{ Gamma prior parm (def: .5)}
    \item{\code{b}}{ Gamma prior parm (def: .1)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
    \item{\code{s_beta}}{ scaling for beta| alpha RW inc cov (def: 2.93/sqrt(nvar))}
    \item{\code{s_alpha}}{ scaling for alpha | beta RW inc cov (def: 2.93)}
    \item{\code{c}}{ fractional likelihood weighting parm (def:2)}
    \item{\code{Vbeta0}}{ starting value for Vbeta (def: I)}
    \item{\code{Delta0}}{ starting value for Delta (def: 0)}
  }
}
\value{
  a list containing: 
  \item{llike}{R/keep vector of values of log-likelihood}
  \item{betadraw}{nreg x nvar x R/keep array of beta draws}
  \item{alphadraw}{R/keep vector of alpha draws}
  \item{acceptrbeta}{acceptance rate of the beta draws}
  \item{acceptralpha}{acceptance rate of the alpha draws}
}
\note{
  The NBD regression encompasses Poisson regression in the sense that as alpha goes to
  infinity the NBD distribution tends to the Poisson.\cr
  For "small" values of alpha, the dependent variable can be extremely variable so that 
  a large number of observations may be required to obtain precise inferences. 

  For ease of interpretation, we recommend demeaning Z variables.
}

\seealso{ \code{\link{rnegbinRw}} }
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 5. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Sridhar Narayanam & Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
##
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

nreg = 100        # Number of cross sectional units
T = 50            # Number of observations per unit
nobs = nreg*T
nvar=2            # Number of X variables
nz=2              # Number of Z variables
              
# Construct the Z matrix
Z = cbind(rep(1,nreg),rnorm(nreg,mean=1,sd=0.125))

Delta = cbind(c(4,2), c(0.1,-1))
alpha = 5
Vbeta = rbind(c(2,1),c(1,2))

# Construct the regdata (containing X)
simnegbindata = NULL
for (i in 1:nreg) {
    betai = as.vector(Z[i,]\%*\%Delta) + chol(Vbeta)\%*\%rnorm(nvar)
    X = cbind(rep(1,T),rnorm(T,mean=2,sd=0.25))
    simnegbindata[[i]] = list(y=simnegbin(X,betai,alpha), X=X,beta=betai)
}

Beta = NULL
for (i in 1:nreg) {Beta=rbind(Beta,matrix(simnegbindata[[i]]$beta,nrow=1))}
    
Data1 = list(regdata=simnegbindata, Z=Z)
Mcmc1 = list(R=R)

out = rhierNegbinRw(Data=Data1, Mcmc=Mcmc1)

cat("Summary of Delta draws",fill=TRUE)
summary(out$Deltadraw,tvalues=as.vector(Delta))
cat("Summary of Vbeta draws",fill=TRUE)
summary(out$Vbetadraw,tvalues=as.vector(Vbeta[upper.tri(Vbeta,diag=TRUE)]))
cat("Summary of alpha draws",fill=TRUE)
summary(out$alpha,tvalues=alpha)

if(0){
## plotting examples
plot(out$betadraw)
plot(out$alpha,tvalues=alpha)
plot(out$Deltadraw,tvalues=as.vector(Delta))
}
}

\keyword{models}
