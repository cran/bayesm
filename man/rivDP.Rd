\name{rivDP}
\alias{rivDP}
\concept{Instrumental Variables}
\concept{Gibbs Sampler}
\concept{Dirichlet Process}
\concept{bayes}
\concept{endogeneity}
\concept{simultaneity}
\concept{MCMC}

\title{ Linear "IV" Model with DP Process Prior for Errors}
\description{
  \code{rivDP} is a Gibbs Sampler for a linear structural equation with an arbitrary number of instruments. 
  \code{rivDP} uses a mixture of normals for the structural and reduced form equation implemented with a 
  Dirichlet Process Prior.
}
\usage{
rivDP(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(z,w,x,y) }
  \item{Prior}{ list(md,Ad,mbg,Abg,lambda,Prioralpha,lambda_hyper) (optional) } 
  \item{Mcmc}{ list(R,keep,nprint,maxuniq,SCALE,gridsize) (R required) }
}
\details{
  Model:\cr
  \eqn{x=z'\delta + e1}. \cr
  \eqn{y=\beta*x + w'\gamma + e2}. \cr
  \eqn{e1,e2} \eqn{\sim}{~} \eqn{N(\theta_{i})}.  \eqn{\theta_{i}} represents \eqn{\mu_{i},\Sigma_{i}}
  
  Note: Error terms have non-zero means.  DO NOT include intercepts in the z or w matrices.  This is different
        from \code{rivGibbs} which requires intercepts to be included explicitly.

  Priors:\cr
  \eqn{\delta} \eqn{\sim}{~} \eqn{N(md,Ad^{-1})}.  \eqn{vec(\beta,\gamma)} \eqn{\sim}{~} \eqn{N(mbg,Abg^{-1})} \cr

  \eqn{\theta_{i}} \eqn{\sim}{~} \eqn{G} \cr

  \eqn{G} \eqn{\sim}{~} \eqn{DP(alpha,G_{0})} \cr
 
  \eqn{G_{0}} is the natural conjugate prior for \eqn{(\mu,\Sigma)}: \cr
  \eqn{\Sigma} \eqn{\sim}{~} \eqn{IW(nu,vI)} and  \eqn{\mu|\Sigma} \eqn{\sim}{~} \eqn{N(0,\Sigma (x) a^{-1})} \cr
  These parameters are collected together in the list \eqn{\lambda}.  It is highly
       recommended that you use the default settings for these hyper-parameters.\cr
  
  \eqn{\lambda(a,nu,v):}\cr
  
      \eqn{a} \eqn{\sim}{~} uniform[alim[1],alimb[2]]\cr
      \eqn{nu} \eqn{\sim}{~} dim(data)-1 + exp(z) \cr
      \eqn{z} \eqn{\sim}{~} uniform[dim(data)-1+nulim[1],nulim[2]]\cr
      \eqn{v} \eqn{\sim}{~} uniform[vlim[1],vlim[2]]
       

  \eqn{alpha} \eqn{\sim}{~} \eqn{(1-(alpha-alpha_{min})/(alpha_{max}-alpha{min}))^{power}} \cr
   where \eqn{alpha_{min}} and \eqn{alpha_{max}} are set using the arguments in the reference
   below.  It is highly recommended that you use the default values for the hyperparameters
   of the prior on alpha

List arguments contain:
  
Data:\cr
  \itemize{
    \item{\code{z}}{ matrix of obs on instruments}
    \item{\code{y}}{ vector of obs on lhs var in structural equation}
    \item{\code{x}}{ "endogenous" var in structural eqn}
    \item{\code{w}}{ matrix of obs on "exogenous" vars in the structural eqn}}
Prior:\cr    
  \itemize{
    \item{\code{md}}{ prior mean of delta (def: 0)}
    \item{\code{Ad}}{ pds prior prec for prior on delta (def: .01I)}
    \item{\code{mbg}}{ prior mean vector for prior on beta,gamma (def: 0)}
    \item{\code{Abg}}{ pds prior prec  for prior on beta,gamma (def: .01I)}}
Prioralpha:\cr
 \itemize{
  \item{\code{Istarmin}}{ expected number of components at lower bound of support of alpha (def: 1)}
  \item{\code{Istarmax}}{ expected number of components at upper bound of support of alpha}
  \item{\code{power}}{ power parameter for alpha prior (def: .8)}
  }
 
lambda_hyper:\cr
  \itemize{
   \item{\code{alim}}{ defines support of a distribution,def:c(.01,10) }
   \item{\code{nulim}}{ defines support of nu distribution, def:c(.01,3)} 
   \item{\code{vlim}}{ defines support of v distribution, def:c(.1,4)} 
  }
  
MCMC:\cr  
  \itemize{
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
    \item{\code{maxuniq}}{ storage constraint on the number of unique components (def: 200)}
    \item{\code{SCALE}}{ scale data (def: TRUE)}
    \item{\code{gridsize}}{ gridsize parm for alpha draws (def: 20)} 
  }

output includes object \code{nmix} of class "bayesm.nmix" which contains draws of predictive distribution of 
  errors (a Bayesian analogue of a density estimate for the error terms).\cr
  nmix:\cr
  \itemize{
    \item{\code{probdraw}}{ not used}
    \item{\code{zdraw}}{ not used}
    \item{\code{compdraw}}{ list R/keep of draws from bivariate predictive for the errors}
  }
  note: in compdraw list, there is only one component per draw
}
\value{
  a list containing:
  \item{deltadraw}{R/keep x dim(delta) array of delta draws}
  \item{betadraw}{R/keep x 1 vector of beta draws}
  \item{gammadraw}{R/keep x dim(gamma) array of gamma draws }
  \item{Istardraw}{R/keep x 1 array of drawsi of the number of unique normal components}
  \item{alphadraw}{R/keep x 1 array of draws of Dirichlet Process tightness parameter}
  \item{nmix}{R/keep x list of draws for predictive distribution of errors}
}
\references{ For further discussion, see "A Semi-Parametric Bayesian Approach to the Instrumental
  Variable Problem," by Conley, Hansen, McCulloch and Rossi, Journal of Econometrics (2008).\cr
}
\seealso{\code{rivGibbs}}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

##
## simulate scaled log-normal errors and run
##
set.seed(66)
k=10
delta=1.5
Sigma=matrix(c(1,.6,.6,1),ncol=2)
N=1000
tbeta=4
set.seed(66)
scalefactor=.6
root=chol(scalefactor*Sigma)
mu=c(1,1)
##
## compute interquartile ranges
##
ninterq=qnorm(.75)-qnorm(.25)
error=matrix(rnorm(100000*2),ncol=2)%*%root
error=t(t(error)+mu)
Err=t(t(exp(error))-exp(mu+.5*scalefactor*diag(Sigma)))
lnNinterq=quantile(Err[,1],prob=.75)-quantile(Err[,1],prob=.25)
##
## simulate data
##
error=matrix(rnorm(N*2),ncol=2)\%*\%root
error=t(t(error)+mu)
Err=t(t(exp(error))-exp(mu+.5*scalefactor*diag(Sigma)))
#
# scale appropriately  
Err[,1]=Err[,1]*ninterq/lnNinterq
Err[,2]=Err[,2]*ninterq/lnNinterq
z=matrix(runif(k*N),ncol=k)
x=z\%*\%(delta*c(rep(1,k)))+Err[,1]
y=x*tbeta+Err[,2]

# set intial values for MCMC
Data = list(); Mcmc=list()
Data$z = z;  Data$x=x; Data$y=y

# start MCMC and keep results
Mcmc$maxuniq=100
Mcmc$R=R
end=Mcmc$R
begin=100

out=rivDP(Data=Data,Mcmc=Mcmc)

cat("Summary of Beta draws",fill=TRUE)
summary(out$betadraw,tvalues=tbeta)

if(0){
## plotting examples
plot(out$betadraw,tvalues=tbeta)
plot(out$nmix)  ## plot "fitted" density of the errors
##

}
}
\keyword{ models }
