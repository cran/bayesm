\name{summary.bayesm.nmix}
\alias{summary.bayesm.nmix}
\title{Summarize Draws of Normal Mixture Components }
\description{
 \code{summary.bayesm.nmix} is an S3 method to display summaries of the distribution implied
 by draws of Normal Mixture Components.  Posterior means and Variance-Covariance matrices are
 displayed.\cr
 \cr
 Note: 1st and 2nd moments may not be very interpretable for mixtures of normals. This summary function
 can take a minute or so.  The current implementation is not efficient.  
}
\usage{
 \method{summary}{bayesm.nmix}(object, names,burnin = trunc(0.1 * nrow(probdraw)), ...)
}
\arguments{
  \item{object}{ an object of class "bayesm.nmix" -- a list of lists of draws}
  \item{names}{ optional character vector of names fo reach dimension of the density}
  \item{burnin}{ number of draws to burn-in (def: .1*nrow(probdraw))}
  \item{...}{ parms to send to summary} 
}
\details{
  an object of class "bayesm.nmix" is a list of three components: 
  \describe{
  \item{probdraw}{ a matrix of R/keep rows by dim of normal mix of mixture prob draws}
  \item{second comp}{ not used}
  \item{compdraw}{ list of list of lists with draws of mixture comp parms}
  }
}
\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}
\seealso{ \code{\link{summary.bayesm.mat}}, \code{\link{summary.bayesm.var}}}
\examples{
##
## not run
#  out=rnmix(Data,Prior,Mcmc)
#  summary(out)
#
}
\keyword{  plot }
