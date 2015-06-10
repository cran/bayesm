\name{rivGibbs}
\alias{rivGibbs}
\concept{Instrumental Variables}
\concept{Gibbs Sampler}
\concept{bayes}
\concept{endogeneity}
\concept{simultaneity}
\concept{MCMC}

\title{ Gibbs Sampler for Linear "IV" Model}
\description{
  \code{rivGibbs} is a Gibbs Sampler for a linear structural equation with an arbitrary number of instruments.
}
\usage{
rivGibbs(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(z,w,x,y) }
  \item{Prior}{ list(md,Ad,mbg,Abg,nu,V) (optional) } 
  \item{Mcmc}{ list(R,keep,nprint) (R required) }
}
\details{
  Model:\cr
  \eqn{x=z'\delta + e1}. \cr
  \eqn{y=\beta*x + w'\gamma + e2}. \cr
  \eqn{e1,e2} \eqn{\sim}{~} \eqn{N(0,\Sigma)}. 
  
  Note: if intercepts are desired in either equation, include vector of ones in z or w

  Priors:\cr
  \eqn{\delta} \eqn{\sim}{~} \eqn{N(md,Ad^{-1})}.  \eqn{vec(\beta,\gamma)} \eqn{\sim}{~} \eqn{N(mbg,Abg^{-1})} \cr
  \eqn{\Sigma} \eqn{\sim}{~} IW(nu,V)

  List arguments contain:
  \itemize{
    \item{\code{z}}{ matrix of obs on instruments}
    \item{\code{y}}{ vector of obs on lhs var in structural equation}
    \item{\code{x}}{ "endogenous" var in structural eqn}
    \item{\code{w}}{ matrix of obs on "exogenous" vars in the structural eqn}
    \item{\code{md}}{ prior mean of delta (def: 0)}
    \item{\code{Ad}}{ pds prior prec for prior on delta (def: .01I)}
    \item{\code{mbg}}{ prior mean vector for prior on beta,gamma (def: 0)}
    \item{\code{Abg}}{ pds prior prec  for prior on beta,gamma (def: .01I)}
    \item{\code{nu}}{ d.f. parm for IW prior on Sigma (def: 5)}
    \item{\code{V}}{ pds location matrix for IW prior on Sigma (def: nuI)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
  }
}
\value{
  a list containing:
  \item{deltadraw}{R/keep x dim(delta) array of delta draws}
  \item{betadraw}{R/keep x 1 vector of beta draws}
  \item{gammadraw}{R/keep x dim(gamma) array of gamma draws }
  \item{Sigmadraw}{R/keep x 4 array of Sigma draws}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 5. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Rob McCulloch and Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(66)
simIV = function(delta,beta,Sigma,n,z,w,gamma) {
eps = matrix(rnorm(2*n),ncol=2) \%*\% chol(Sigma)
x = z \%*\% delta + eps[,1]; y = beta*x +  eps[,2] + w\%*\%gamma
list(x=as.vector(x),y=as.vector(y)) }
n = 200 ; p=1 # number of instruments
z = cbind(rep(1,n),matrix(runif(n*p),ncol=p))
w = matrix(1,n,1)
rho=.8
Sigma = matrix(c(1,rho,rho,1),ncol=2)
delta = c(1,4); beta = .5; gamma = c(1)
simiv = simIV(delta,beta,Sigma,n,z,w,gamma)

Mcmc1=list();  Data1 = list()
Data1$z = z; Data1$w=w; Data1$x=simiv$x; Data1$y=simiv$y
Mcmc1$R = R
Mcmc1$keep=1
out=rivGibbs(Data=Data1,Mcmc=Mcmc1)

cat("Summary of Beta draws",fill=TRUE)
summary(out$betadraw,tvalues=beta)
cat("Summary of Sigma draws",fill=TRUE)
summary(out$Sigmadraw,tvalues=as.vector(Sigma[upper.tri(Sigma,diag=TRUE)]))

if(0){
## plotting examples
plot(out$betadraw)
}
}
\keyword{ models }
