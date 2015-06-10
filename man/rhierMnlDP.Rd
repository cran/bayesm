\name{rhierMnlDP}
\alias{rhierMnlDP}
\concept{bayes}
\concept{MCMC}
\concept{Multinomial Logit}
\concept{normal mixture}
\concept{Dirichlet Process Prior}
\concept{heterogeneity}
\concept{hierarchical models}

\title{ MCMC Algorithm for Hierarchical Multinomial Logit with Dirichlet Process Prior Heterogeneity}
\description{
  \code{rhierMnlDP} is a MCMC algorithm for a hierarchical multinomial logit with a Dirichlet Process Prior for the distribution of heteorogeneity.  A base normal model is used so that the DP can be interpreted as allowing for a mixture of normals with as many components as there are panel units.  This is a hybrid Gibbs Sampler with a RW Metropolis step for the MNL 
  coefficients for each panel unit.  This procedure can be interpreted as a Bayesian semi-parameteric method in the sense that the DP prior can accomodate heterogeniety of an unknown form.
}
\usage{
rhierMnlDP(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(p,lgtdata,Z) ( Z is optional) }
  \item{Prior}{ list(deltabar,Ad,Prioralpha,lambda_hyper) (all are optional)}
  \item{Mcmc}{ list(s,w,R,keep,nprint) (R required)}
}
\details{
  Model: \cr
  \eqn{y_i} \eqn{\sim}{~} \eqn{MNL(X_i,\beta_i)}.  i=1,\ldots, length(lgtdata). \eqn{\theta_i} is nvar x 1.

  \eqn{\beta_i}= Z\eqn{\Delta}[i,] + \eqn{u_i}. \cr
  Note:  Z\eqn{\Delta} is the matrix Z * \eqn{\Delta}; [i,] refers to ith row of this product.\cr
  Delta is an nz x nvar array. 

  \eqn{\beta_i} \eqn{\sim}{~} \eqn{N(\mu_i,\Sigma_i)}. \cr

  Priors: \cr
        \eqn{\theta_i=(\mu_i,\Sigma_i)} \eqn{\sim}{~} \eqn{DP(G_0(\lambda),alpha)}\cr
        \eqn{G_0(\lambda):}\cr
        \eqn{\mu_i | \Sigma_i} \eqn{\sim}{~} \eqn{N(0,\Sigma_i (x) a^{-1})}\cr
        \eqn{\Sigma_i} \eqn{\sim}{~} \eqn{IW(nu,nu*v*I)}\cr
        \eqn{delta= vec(\Delta)} \eqn{\sim}{~} \eqn{N(deltabar,A_d^{-1})}\cr
        
        \eqn{\lambda(a,nu,v):}\cr
        \eqn{a} \eqn{\sim}{~} uniform[alim[1],alimb[2]]\cr
        \eqn{nu} \eqn{\sim}{~}  dim(data)-1 + exp(z) \cr
        \eqn{z} \eqn{\sim}{~}  uniform[dim(data)-1+nulim[1],nulim[2]]\cr
        \eqn{v} \eqn{\sim}{~} uniform[vlim[1],vlim[2]]
       
        \eqn{alpha} \eqn{\sim}{~} \eqn{(1-(alpha-alphamin)/(alphamax-alphamin))^{power}} \cr
        alpha = alphamin then expected number of components = Istarmin \cr
        alpha = alphamax then expected number of components = Istarmax \cr


  Lists contain: \cr

Data:\cr
  \itemize{
    \item{\code{p}}{ p is number of choice alternatives}
    \item{\code{lgtdata}}{list of lists with each cross-section unit MNL data}
    \item{\code{lgtdata[[i]]$y}}{ \eqn{n_i} vector of multinomial outcomes (1,\ldots,m)}
    \item{\code{lgtdata[[i]]$X}}{ \eqn{n_i} by nvar design matrix for ith unit}
  }
Prior: \cr
    \itemize{
      \item{\code{deltabar}}{nz*nvar vector of prior means (def: 0)}
      \item{\code{Ad}}{ prior prec matrix for vec(D) (def: .01I)}
   }
Prioralpha:\cr
 \itemize{
  \item{\code{Istarmin}}{expected number of components at lower bound of support of alpha def(1)}
  \item{\code{Istarmax}}{expected number of components at upper bound of support of alpha (def: min(50,.1*nlgt))}
  \item{\code{power}}{power parameter for alpha prior (def: .8)}
  }
 
lambda_hyper:\cr
  \itemize{
   \item{\code{alim}}{defines support of a distribution (def: (.01,2))}
   \item{\code{nulim}}{defines support of nu distribution (def: (.01,3))} 
   \item{\code{vlim}}{defines support of v distribution (def: (.1,4))} 
  }

Mcmc:\cr
 \itemize{
   \item{\code{R}}{ number of mcmc draws}
   \item{\code{keep}}{ thinning parm, keep every keepth draw}
   \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
   \item{\code{maxuniq}}{ storage constraint on the number of unique components}
   \item{\code{gridsize}}{ number of discrete points for hyperparameter priors,def: 20}
  }

}
\value{
  a list containing:
  \item{Deltadraw}{R/keep  x nz*nvar matrix of draws of Delta, first row is initial value}
  \item{betadraw}{ nlgt x nvar x R/keep array of draws of betas}
  \item{nmix}{ list of 3 components, probdraw, NULL, compdraw }
  \item{adraw}{R/keep draws of hyperparm a}
  \item{vdraw}{R/keep draws of hyperparm v}
  \item{nudraw}{R/keep draws of hyperparm nu}
  \item{Istardraw}{R/keep draws of number of unique components}
  \item{alphadraw}{R/keep draws of number of DP tightness parameter}
  \item{loglike}{R/keep draws of log-likelihood}
}
\note{

  As is well known, Bayesian density estimation involves computing the predictive distribution of a "new" unit parameter,
  \eqn{\theta_{n+1}} (here "n"=nlgt). This is done by averaging the normal base distribution over draws from the distribution of \eqn{\theta_{n+1}} given \eqn{\theta_1}, ..., \eqn{\theta_n},alpha,lambda,Data.
  To facilitate this, we store those draws from the predictive distribution of \eqn{\theta_{n+1}} in a list structure compatible  with other \code{bayesm} routines that implement a finite mixture of normals.

  More on nmix list:\cr 
  contains the draws from the predictive distribution of a "new" observations parameters.  These are simply the parameters of one normal distribution.  We enforce compatibility with a mixture of k components in order to utilize generic summary 
  plotting functions.  

  Therefore,\code{probdraw} is a vector of ones.  \code{zdraw} (indicator draws) is omitted as it is not necessary for density estimation. \code{compdraw} contains the draws of the \eqn{\theta_{n+1}} as a list of list of lists.

  More on \code{compdraw} component of return value list:
  \itemize{
  \item{compdraw[[i]]}{ith draw of components for mixtures}
  \item{compdraw[[i]][[1]]}{ith draw of the thetanp1}
  \item{compdraw[[i]][[1]][[1]]}{ith draw of mean vector}
  \item{compdraw[[i]][[1]][[2]]}{ith draw of parm (rooti)}
  }

  We parameterize the prior on \eqn{\Sigma_i} such that \eqn{mode(\Sigma)= nu/(nu+2) vI}.
    The support of nu enforces a non-degenerate IW density; \eqn{nulim[1] > 0}.

    The default choices of alim,nulim, and vlim determine the location and approximate size of candidate
    "atoms" or possible normal components. The defaults are sensible given a reasonable scaling of the X variables.
    You want to insure that alim is set for a wide enough range of values (remember a is a precision
    parameter) and the v is big enough to propose Sigma matrices wide enough to cover the data range.  

    A careful analyst should look at the posterior distribution of a, nu, v to make sure that the support is
    set correctly in alim, nulim, vlim.  In other words, if we see the posterior bunched up at one end of these
    support ranges, we should widen the range and rerun.  

   If you want to force the procedure to use many small atoms, then set nulim to consider only large values and 
   set vlim to consider only small scaling constants.  Set alphamax to a large number.  This will create a very
   "lumpy" density estimate somewhat like the classical Kernel density estimates. Of course, this is not advised 
   if you have a prior belief that densities are relatively smooth.

  Note: Z should \strong{not} include an intercept and is centered for ease of interpretation. The mean of each of the \code{nlgt} \eqn{\beta} s is the mean of the normal mixture.  Use \code{summary()} to compute this mean from the \code{compdraw} output.\cr
  
  Large R values may be required (>20,000).

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
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=20000} else {R=10}

set.seed(66)
p=3                                # num of choice alterns
ncoef=3  
nlgt=300                           # num of cross sectional units
nz=2
Z=matrix(runif(nz*nlgt),ncol=nz)
Z=t(t(Z)-apply(Z,2,mean))          # demean Z
ncomp=3                                # no of mixture components
Delta=matrix(c(1,0,1,0,1,2),ncol=2)
comps=NULL
comps[[1]]=list(mu=c(0,-1,-2),rooti=diag(rep(2,3)))
comps[[2]]=list(mu=c(0,-1,-2)*2,rooti=diag(rep(2,3)))
comps[[3]]=list(mu=c(0,-1,-2)*4,rooti=diag(rep(2,3)))
pvec=c(.4,.2,.4)

simmnlwX= function(n,X,beta) {
  ##  simulate from MNL model conditional on X matrix
  k=length(beta)
  Xbeta=X\%*\%beta
  j=nrow(Xbeta)/n
  Xbeta=matrix(Xbeta,byrow=TRUE,ncol=j)
  Prob=exp(Xbeta)
  iota=c(rep(1,j))
  denom=Prob\%*\%iota
  Prob=Prob/as.vector(denom)
  y=vector("double",n)
  ind=1:j
  for (i in 1:n) 
      {yvec=rmultinom(1,1,Prob[i,]); y[i]=ind\%*\%yvec}
  return(list(y=y,X=X,beta=beta,prob=Prob))
}

## simulate data with a mixture of 3 normals
simlgtdata=NULL
ni=rep(50,300)
for (i in 1:nlgt) 
{  betai=Delta\%*\%Z[i,]+as.vector(rmixture(1,pvec,comps)$x)
   Xa=matrix(runif(ni[i]*p,min=-1.5,max=0),ncol=p)
   X=createX(p,na=1,nd=NULL,Xa=Xa,Xd=NULL,base=1)
   outa=simmnlwX(ni[i],X,betai)
   simlgtdata[[i]]=list(y=outa$y,X=X,beta=betai)
}

## plot betas
if(1){
## set if(1) above to produce plots
bmat=matrix(0,nlgt,ncoef)
for(i in 1:nlgt) {bmat[i,]=simlgtdata[[i]]$beta}
par(mfrow=c(ncoef,1))
for(i in 1:ncoef) hist(bmat[,i],breaks=30,col="magenta")
}

##   set Data and Mcmc lists
keep=5
Mcmc1=list(R=R,keep=keep)
Data1=list(p=p,lgtdata=simlgtdata,Z=Z)

out=rhierMnlDP(Data=Data1,Mcmc=Mcmc1)

cat("Summary of Delta draws",fill=TRUE)
summary(out$Deltadraw,tvalues=as.vector(Delta))

if(0) {
## plotting examples
plot(out$betadraw)
plot(out$nmix)
}

}

\keyword{models}
