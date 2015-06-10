\name{plot.bayesm.nmix}
\alias{plot.bayesm.nmix}
\concept{MCMC}
\concept{S3 method}
\concept{plot}
\title{Plot Method for MCMC Draws of Normal Mixtures}
\description{
   \code{plot.bayesm.nmix} is an S3 method to plot aspects of the fitted density from a list
    of MCMC draws of normal mixture components. Plots of marginal univariate and bivariate densities
    are produced.
}
\usage{
\method{plot}{bayesm.nmix}(x,names,burnin,Grid,bi.sel,nstd,marg,Data,ngrid,ndraw, ...)
}
\arguments{
  \item{x}{ An object of  S3 class bayesm.nmix }
  \item{names}{optional character vector of names for each of the dimensions}
  \item{burnin}{number of draws to discard for burn-in (def: .1*nrow(X))}
  \item{Grid}{matrix of grid points for densities, def: mean +/- nstd std deviations (if Data no supplied), 
               range of Data if supplied)}
  \item{bi.sel}{list of vectors, each giving pairs for bivariate distributions (def: list(c(1,2)))}
  \item{nstd}{number of standard deviations for default Grid (def: 2)}
  \item{marg}{logical, if TRUE display marginals (def: TRUE)}
  \item{Data}{matrix of data points, used to paint histograms on marginals and for grid  }
  \item{ngrid}{number of grid points for density estimates (def: 50)}
  \item{ndraw}{number of draws to average Mcmc estimates over (def: 200)}
  \item{...}{ standard graphics parameters }
}
\details{
  Typically, \code{plot.bayesm.nmix} will be invoked by a call to the generic plot function as in
  \code{plot(object)} where object is of class bayesm.nmix.  These objects are lists of three components. 
  The first component is an array of draws of mixture component probabilties.  The second component
  is not used. The third is a lists of lists of lists with draws of each of the normal components.
  \cr
  \cr
  \code{plot.bayesm.nmix} can also be used as a standard function, as in \code{plot.bayesm.nmix(list)}.
}
\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\seealso{ \code{\link{rnmixGibbs}}, \code{\link{rhierMnlRwMixture}}, \code{\link{rhierLinearMixture}},
           \code{\link{rDPGibbs}}}
\examples{
##
## not run
#  out=rnmixGibbs(Data,Prior,Mcmc)
#  plot(out,bi.sel=list(c(1,2),c(3,4),c(1,3)))
#        # plot bivariate distributions for dimension 1,2; 3,4; and 1,3
#

}
\keyword{  hplot }
