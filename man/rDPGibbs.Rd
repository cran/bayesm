\name{rDPGibbs}
\alias{rDPGibbs}
\concept{bayes}
\concept{MCMC}
\concept{normal mixtures}
\concept{Dirichlet Process}
\concept{Gibbs Sampling}
\title{ Density Estimation with Dirichlet Process Prior and Normal Base }
\description{
    \code{rDPGibbs} implements a Gibbs Sampler to draw from the posterior for a normal mixture problem
    with a Dirichlet Process prior.  A natural conjugate base prior is used along with priors on the hyper 
    parameters of this distribution. One interpretation
    of this model is as a normal mixture with a random number of components that can grow with the sample size. 
}

\usage{
rDPGibbs(Prior, Data, Mcmc)
}

\arguments{
  \item{Prior}{ list(Prioralpha,lambda_hyper) }
  \item{Data}{ list(y) }
  \item{Mcmc}{ list(R,keep,nprint,maxuniq,SCALE,gridsize) }
}

\details{

Model: \cr
        \eqn{y_i} \eqn{\sim}{~} \eqn{N(\mu_i,\Sigma_i)}. \cr

Priors:\cr
        \eqn{\theta_i=(\mu_i,\Sigma_i)} \eqn{\sim}{~} \eqn{DP(G_0(\lambda),alpha)}\cr
        \eqn{G_0(\lambda):}\cr
        \eqn{\mu_i | \Sigma_i} \eqn{\sim}{~} \eqn{N(0,\Sigma_i (x) a^{-1})}\cr
        \eqn{\Sigma_i} \eqn{\sim}{~} \eqn{IW(nu,nu*v*I)}
        
        \eqn{\lambda(a,nu,v):}\cr
        \eqn{a} \eqn{\sim}{~} uniform on grid[alim[1],alimb[2]]\cr
        \eqn{nu} \eqn{\sim}{~} uniform on grid[dim(data)-1 + exp(nulim[1]),dim(data)-1 +exp(nulim[2])]\cr
        \eqn{v} \eqn{\sim}{~} uniform on grid[vlim[1],vlim[2]]
       
        \eqn{alpha} \eqn{\sim}{~} \eqn{(1-(\alpha-alphamin)/(alphamax-alphamin))^{power}} \cr
        \eqn{alpha}= alphamin then expected number of components = Istarmin \cr
        \eqn{alpha}= alphamax then expected number of components = Istarmax \cr

List arguments contain:

Data:\cr
  \itemize{
    \item{\code{y}}{N x k matrix of observations on k dimensional data}
  }

Prioralpha:\cr
 \itemize{
  \item{\code{Istarmin}}{ expected number of components at lower bound of support of alpha (def: 1)}
  \item{\code{Istarmax}}{ expected number of components at upper bound of support of alpha}
  \item{\code{power}}{ power parameter for alpha prior (def: .8)}
  }
 
lambda_hyper:\cr
  \itemize{
   \item{\code{alim}}{ defines support of a distribution (def: (.01,10))}
   \item{\code{nulim}}{ defines support of nu distribution (def: (.01,3))} 
   \item{\code{vlim}}{ defines support of v distribution (def: (.1,4))} 
  }
Mcmc:\cr
 \itemize{
   \item{\code{R}}{ number of mcmc draws}
   \item{\code{keep}}{ thinning parm, keep every keepth draw}
   \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}
   \item{\code{maxuniq}}{ storage constraint on the number of unique components (def: 200)}
   \item{\code{SCALE}}{ should data be scaled by mean,std deviation before posterior draws, (def: TRUE)}
   \item{\code{gridsize}}{ number of discrete points for hyperparameter priors,def: 20}
  }

the basic output are draws from the predictive distribution of the data in the object, \code{nmix}. 
The average of these draws is the Bayesian analogue of a density estimate.

nmix:\cr
  \itemize{
   \item{\code{probdraw}}{ R/keep x 1 matrix of 1s}
   \item{\code{zdraw}}{ R/keep x N matrix of draws of indicators of which component each obs is assigned to}
   \item{\code{compdraw}}{ R/keep list of draws of normals}
  }
  Output of the components is in the form of a list of lists. \cr
  compdraw[[i]] is ith draw -- list of lists. \cr
  compdraw[[i]][[1]] is list of parms for a draw from predictive. \cr
  compdraw[[i]][1]][[1]] is the mean vector. compdraw[[i]][[1]][[2]] is the inverse of Cholesky root.
  \eqn{\Sigma} = t(R)\%*\%R, \eqn{R^{-1}} = compdraw[[i]][[1]][[2]].
}


\note{
    we parameterize the prior on \eqn{\Sigma_i} such that \eqn{mode(\Sigma)= nu/(nu+2) vI}.
    The support of nu enforces valid IW density; \eqn{nulim[1] > 0}

    We use the structure for \code{nmix} that is compatible with the \code{bayesm} routines for finite mixtures of normals.
    This allows us to use the same summary and plotting methods.  

    The default choices of alim,nulim, and vlim determine the location and approximate size of candidate
    "atoms" or possible normal components. The defaults are sensible given that we scale the data.  Without 
    scaling, you want to insure that alim is set for a wide enough range of values (remember a is a precision
    parameter) and the v is big enough to propose Sigma matrices wide enough to cover the data range.  

    A careful analyst should look at the posterior distribution of a, nu, v to make sure that the support is
    set correctly in alim, nulim, vlim.  In other words, if we see the posterior bunched up at one end of these
    support ranges, we should widen the range and rerun.  

   If you want to force the procedure to use many small atoms, then set nulim to consider only large values and 
   set vlim to consider only small scaling constants.  Set Istarmax to a large number.  This will create a very
   "lumpy" density estimate somewhat like the classical Kernel density estimates. Of course, this is not advised 
   if you have a prior belief that densities are relatively smooth.

}


\value{
 \item{nmix}{ a list containing: probdraw,zdraw,compdraw}
 \item{alphadraw}{ vector of draws of DP process tightness parameter}
 \item{nudraw}{ vector of draws of base prior hyperparameter}
 \item{adraw}{ vector of draws of base prior hyperparameter}
 \item{vdraw}{ vector of draws of base prior hyperparameter}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\seealso{ \code{\link{rnmixGibbs}},\code{\link{rmixture}}, \code{\link{rmixGibbs}} ,
          \code{\link{eMixMargDen}}, \code{\link{momMix}}, \code{\link{mixDen}}, \code{\link{mixDenBi}}}

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

## simulate univariate data from Chi-Sq

set.seed(66)
N=200
chisqdf=8; y1=as.matrix(rchisq(N,df=chisqdf))

## set arguments for rDPGibbs

Data1=list(y=y1)
Prioralpha=list(Istarmin=1,Istarmax=10,power=.8)
Prior1=list(Prioralpha=Prioralpha)

Mcmc=list(R=R,keep=1,maxuniq=200)

out1=rDPGibbs(Prior=Prior1,Data=Data1,Mcmc)

if(0){
## plotting examples
rgi=c(0,20); grid=matrix(seq(from=rgi[1],to=rgi[2],length.out=50),ncol=1)
deltax=(rgi[2]-rgi[1])/nrow(grid)
plot(out1$nmix,Grid=grid,Data=y1)
## plot true density with historgram
plot(range(grid[,1]),1.5*range(dchisq(grid[,1],df=chisqdf)),
  type="n",xlab=paste("Chisq ; ",N," obs",sep=""), ylab="")
hist(y1,xlim=rgi,freq=FALSE,col="yellow",breaks=20,add=TRUE)
lines(grid[,1],dchisq(grid[,1],df=chisqdf)/
  (sum(dchisq(grid[,1],df=chisqdf))*deltax),col="blue",lwd=2)
}


## simulate bivariate data from the  "Banana" distribution (Meng and Barnard) 
banana=function(A,B,C1,C2,N,keep=10,init=10)
{ R=init*keep+N*keep
   x1=x2=0
   bimat=matrix(double(2*N),ncol=2)
  for (r in 1:R)
  { x1=rnorm(1,mean=(B*x2+C1)/(A*(x2^2)+1),sd=sqrt(1/(A*(x2^2)+1)))
  x2=rnorm(1,mean=(B*x2+C2)/(A*(x1^2)+1),sd=sqrt(1/(A*(x1^2)+1)))
  if (r>init*keep && r\%\%keep==0) {mkeep=r/keep; bimat[mkeep-init,]=c(x1,x2)} }
return(bimat)
}


set.seed(66)
nvar2=2
A=0.5; B=0; C1=C2=3
y2=banana(A=A,B=B,C1=C1,C2=C2,1000)

Data2=list(y=y2)
Prioralpha=list(Istarmin=1,Istarmax=10,power=.8)
Prior2=list(Prioralpha=Prioralpha)
Mcmc=list(R=R,keep=1,maxuniq=200)

out2=rDPGibbs(Prior=Prior2,Data=Data2,Mcmc)


if(0){
## plotting examples

rx1=range(y2[,1]); rx2=range(y2[,2])
x1=seq(from=rx1[1],to=rx1[2],length.out=50)
x2=seq(from=rx2[1],to=rx2[2],length.out=50)
grid=cbind(x1,x2)

plot(out2$nmix,Grid=grid,Data=y2)

## plot true bivariate density
tden=matrix(double(50*50),ncol=50)
for (i in 1:50){ for (j in 1:50) 
      {tden[i,j]=exp(-0.5*(A*(x1[i]^2)*(x2[j]^2)+
      (x1[i]^2)+(x2[j]^2)-2*B*x1[i]*x2[j]-2*C1*x1[i]-2*C2*x2[j]))}
}
tden=tden/sum(tden)
image(x1,x2,tden,col=terrain.colors(100),xlab="",ylab="")
contour(x1,x2,tden,add=TRUE,drawlabels=FALSE)
title("True Density")
}
}
\keyword{ multivariate }
