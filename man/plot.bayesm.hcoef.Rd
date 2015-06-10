\name{plot.bayesm.hcoef}
\alias{plot.bayesm.hcoef}
\concept{MCMC}
\concept{S3 method}
\concept{plot}
\concept{hierarchical model}
\title{Plot Method for Hierarchical Model Coefs }
\description{
   \code{plot.bayesm.hcoef} is an S3 method to plot 3 dim arrays of hierarchical coefficients. 
   Arrays are of class bayesm.hcoef with dimensions: cross-sectional unit x coef x MCMC draw. 
}
\usage{
\method{plot}{bayesm.hcoef}(x,names,burnin,...)
}
\arguments{
  \item{x}{ An object of S3 class, bayesm.hcoef }
  \item{names}{ a list of names for the variables in the hierarchical model}
  \item{burnin}{ no draws to burnin (def: .1*R)}
  \item{...}{ standard graphics parameters }
}
\details{
  Typically, \code{plot.bayesm.hcoef} will be invoked by a call to the generic plot function as in
  \code{plot(object)} where object is of class bayesm.hcoef.  All of the \code{bayesm} hierarchical routines
  return draws of hierarchical coefficients in this class (see example below).  One can also simply invoke 
  \code{plot.bayesm.hcoef} on any valid 3-dim array as in \code{plot.bayesm.hcoef(betadraws)}
  \cr
  \cr
  \code{plot.bayesm.hcoef} is also exported for use as a standard function, as in
       \code{plot.bayesm.hcoef(array)}.
}
\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\seealso{ \code{\link{rhierMnlRwMixture}},\code{\link{rhierLinearModel}},
          \code{\link{rhierLinearMixture}},\code{\link{rhierNegbinRw}} }
\examples{
##
## not run
#  out=rhierLinearModel(Data,Prior,Mcmc)
#  plot(out$betadraws)
#
}
\keyword{ hplot }
