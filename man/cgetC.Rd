\name{cgetC}
\alias{cgetC}

\title{Obtain A List of Cut-offs for Scale Usage Problems}

\description{
\code{cgetC} obtains a list of censoring points, or cut-offs, used in the ordinal multivariate probit model of Rossi et al (2001). This approach uses a quadratic parameterization of the cut-offs. The model is useful for modeling correlated ordinal data on a scale from \eqn{1} to \eqn{k} with different scale usage patterns. 
}

\usage{cgetC(e, k)}

\arguments{
  \item{e}{ quadratic parameter (\eqn{0 <} \code{e} \eqn{< 1}) }
  \item{k}{ items are on a scale from \eqn{1,\ldots, k} }
}

\section{Warning}{This is a utility function which implements \strong{no} error-checking.}

\value{A vector of \eqn{k+1} cut-offs.}

\author{Rob McCulloch and Peter Rossi, Anderson School, UCLA. \email{perossichi@gmail.com}.}

\references{Rossi et al (2001), \dQuote{Overcoming Scale Usage Heterogeneity,} \emph{JASA} 96, 20--31.}

\seealso{ \code{\link{rscaleUsage}} }

\examples{
cgetC(0.1, 10)
}

\keyword{ utilities }
