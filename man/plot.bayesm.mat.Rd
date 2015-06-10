\name{plot.bayesm.mat}
\alias{plot.bayesm.mat}
\concept{MCMC}
\concept{S3 method}
\concept{plot}
\title{Plot Method for Arrays of MCMC Draws}
\description{
   \code{plot.bayesm.mat} is an S3 method to plot arrays of MCMC draws. The columns in the
   array correspond to parameters and the rows to MCMC draws.
}
\usage{
\method{plot}{bayesm.mat}(x,names,burnin,tvalues,TRACEPLOT,DEN,INT,CHECK_NDRAWS, ...)
}
\arguments{
  \item{x}{ An object of either S3 class, bayesm.mat, or S3 class, mcmc }
  \item{names}{optional character vector of names for coefficients}
  \item{burnin}{number of draws to discard for burn-in (def: .1*nrow(X))}
  \item{tvalues}{vector of true values}
  \item{TRACEPLOT}{ logical, TRUE provide sequence plots of draws and acfs (def: TRUE)}
  \item{DEN}{ logical, TRUE use density scale on histograms (def: TRUE)}
  \item{INT}{ logical, TRUE put various intervals and points on graph (def: TRUE)}
  \item{CHECK_NDRAWS}{ logical, TRUE check that there are at least 100 draws (def: TRUE)}
  \item{...}{ standard graphics parameters }
}
\details{
  Typically, \code{plot.bayesm.mat} will be invoked by a call to the generic plot function as in
  \code{plot(object)} where object is of class bayesm.mat.  All of the \code{bayesm} MCMC routines
  return draws in this class (see example below).  One can also simply invoke 
  \code{plot.bayesm.mat} on any valid 2-dim array as in \code{plot.bayesm.mat(betadraws)}. \cr
  \cr
  \code{plot.bayesm.mat} paints (by default) on the histogram: \cr
  \cr
  green "[]" delimiting 95\% Bayesian Credibility Interval \cr
  yellow "()" showing +/- 2 numerical standard errors \cr
  red "|" showing posterior mean 
  \cr
  \cr
  \code{plot.bayesm.mat} is also exported for use as a standard function, as in
      \code{plot.bayesm.mat(matrix)}
}
\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\examples{
##
## not run
#  out=runiregGibbs(Data,Prior,Mcmc)
#  plot(out$betadraw)
#
}
\keyword{ hplot  }

