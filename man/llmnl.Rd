\name{llmnl}
\alias{llmnl}

\title{ Evaluate Log Likelihood for Multinomial Logit Model  }
\description{
  \code{llmnl} evaluates log-likelihood for the multinomial logit model.
}
\usage{
llmnl(y, X, beta)
}

\arguments{
  \item{y}{ n x 1 vector of obs on y (1,\ldots, p) }
  \item{X}{ n*p x k Design matrix (use \code{createX} to make) }
  \item{beta}{ k x 1 coefficient vector }
}
\details{
  Let \eqn{mu_i=X_i}\%*\%\eqn{\beta}, then \eqn{Pr(y=j) = e^{mu_j}/\sum_i{mu_i}}.\cr
  \eqn{X_i} is the submatrix of X corresponding to the
  ith observation.  X has n*p rows.  
  
  Use \code{\link{createX}} to create X. 
}
\value{
  value of log-likelihood (sum of log prob of observed multinomial outcomes).
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{createX}}, \code{\link{rmnlIndepMetrop}} }

\examples{
##
\dontrun{ll=llmnl(y,X,beta)}
}

\keyword{ models }
