\name{clusterMix}
\alias{clusterMix}
\concept{normal mixture}
\concept{clustering}
\title{ Cluster Observations Based on Indicator MCMC Draws }
\description{
  \code{clusterMix} uses MCMC draws of indicator variables from a normal
  component mixture model to cluster observations based on a similarity matrix.
}
\usage{
clusterMix(zdraw, cutoff = 0.9, SILENT = FALSE, nprint = BayesmConstant.nprint)
}
\arguments{
  \item{zdraw}{ R x nobs array of draws of indicators }
  \item{cutoff}{ cutoff probability for similarity  (def: .9)}
  \item{SILENT}{ logical flag for silent operation (def: FALSE) }
  \item{nprint}{ print every nprint'th draw (def: 100) }
}
\details{

   Define a similarity matrix, Sim, Sim[i,j]=1 if observations i and j are in same component.
   Compute the posterior mean of Sim over indicator draws.

   Clustering is achieved by two means:

   Method A:
   Find the indicator draw whose similarity matrix minimizes, loss(E[Sim]-Sim(z)),  
   where loss is absolute deviation.

   Method B:
   Define a Similarity matrix by setting any element of E[Sim] = 1 if E[Sim] > cutoff.
   Compute the clustering scheme associated with this "windsorized" Similarity matrix.
}
\value{
  \item{clustera}{indicator function for clustering based on method A above}
  \item{clusterb}{indicator function for clustering based on method B above}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch Chapter 3. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA, <perossichi@gmail.com>.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}
\seealso{ \code{\link{rnmixGibbs}}  }

\keyword{ models }
\keyword{ multivariate }
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) 
{
## simulate data from mixture of normals
n=500
pvec=c(.5,.5)
mu1=c(2,2)
mu2=c(-2,-2)
Sigma1=matrix(c(1,.5,.5,1),ncol=2)
Sigma2=matrix(c(1,.5,.5,1),ncol=2)
comps=NULL
comps[[1]]=list(mu1,backsolve(chol(Sigma1),diag(2)))
comps[[2]]=list(mu2,backsolve(chol(Sigma2),diag(2)))
dm=rmixture(n,pvec,comps)
## run MCMC on normal mixture
R=2000
Data=list(y=dm$x)
ncomp=2
Prior=list(ncomp=ncomp,a=c(rep(100,ncomp)))
Mcmc=list(R=R,keep=1)
out=rnmixGibbs(Data=Data,Prior=Prior,Mcmc=Mcmc)
begin=500
end=R
## find clusters
outclusterMix=clusterMix(out$nmix$zdraw[begin:end,])
##
## check on clustering versus "truth"
##  note: there could be switched labels
##
table(outclusterMix$clustera,dm$z)
table(outclusterMix$clusterb,dm$z)
}
##
}
