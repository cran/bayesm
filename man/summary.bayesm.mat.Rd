\name{summary.bayesm.mat}
\alias{summary.bayesm.mat}
\title{Summarize Mcmc Parameter Draws }
\description{
 \code{summary.bayesm.mat} is an S3 method to summarize marginal distributions given an
  array of draws
}
\usage{
\method{summary}{bayesm.mat}(object, names, burnin = trunc(0.1 * nrow(X)), 
  tvalues, QUANTILES = TRUE, TRAILER = TRUE,...)
}
\arguments{
  \item{object}{ \code{object} (hereafter \code{X}) is an array of draws, usually an object of class "bayesm.mat" }
  \item{names}{ optional character vector of names for the columns of \code{X}}
  \item{burnin}{ number of draws to burn-in (def: .1*nrow(X))}
  \item{tvalues}{ optional vector of "true" values for use in simulation examples }
  \item{QUANTILES}{ logical for should quantiles be displayed (def: TRUE)}
  \item{TRAILER}{ logical for should a trailer be displayed (def: TRUE)}
  \item{...}{ optional arguments for generic function }
}
\details{
  Typically, \code{summary.bayesm.nmix} will be invoked by a call to the generic summary function as in
  \code{summary(object)} where object is of class bayesm.mat. Mean, Std Dev, Numerical Standard error (of 
  estimate of posterior mean), relative numerical efficiency (see \code{numEff}) and effective sample
  size are displayed.  If QUANTILES=TRUE, quantiles of marginal distirbutions in the columns of X are displayed.  \cr
  \cr
  \code{summary.bayesm.mat} is also exported for direct use as a standard function, as in 
        \code{summary.bayesm.mat(matrix)}.
   \cr
   \code{summary.bayesm.mat(matrix)} returns (invisibly) the array of the various summary statistics for further use.  To 
   assess this array use\code{stats=summary(Drawmat)}.
}
\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\seealso{ \code{\link{summary.bayesm.var}}, \code{\link{summary.bayesm.nmix}}}
\examples{
##
## not run
#  out=rmnpGibbs(Data,Prior,Mcmc)
#  summary(out$betadraw)
#

}
\keyword{ univar }

