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
\method{plot}{bayesm.nmix}(x,names,burnin,Grid,bi.sel,nstd, ...)
}
\arguments{
  \item{x}{ An object of either S3 class, bayesm.nmix, or S3 class, mcmc }
  \item{names}{optional character vector of names for each of the dimensions}
  \item{burnin}{number of draws to discard for burn-in, def: .1*nrow(X)}
  \item{Grid}{matrix of grid points for densities, def: mean +/- nstd std deviations}
  \item{bi.sel}{list of vectors, each giving pairs for bivariate distributions, def: list(c(1,2))}
  \item{nstd}{number of standard deviations for default Grid, def: 2}
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
\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}
\seealso{ \code{\link{rnmixGibbs}}, \code{\link{rhierMnlRwMixture}}, \code{\link{rhierLinearMixture}}}
\examples{
##
## not run
#  out=rnmix(Data,Prior,Mcmc)
#  plot(out,bi.sel=list(c(1,2),c(3,4),c(1,3)))
#        plot bivariate distributions for dimension 1,2; 3,4; and 1,3
#

}
\keyword{  hplot }
