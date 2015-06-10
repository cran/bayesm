\name{llmnl}
\alias{llmnl}
\concept{multinomial logit}
\concept{likelihood}

\title{ Evaluate Log Likelihood for Multinomial Logit Model  }
\description{
  \code{llmnl} evaluates log-likelihood for the multinomial logit model.
}
\usage{
llmnl(beta,y, X)
}

\arguments{
  \item{beta}{ k x 1 coefficient vector }
  \item{y}{ n x 1 vector of obs on y (1,\ldots, p) }
  \item{X}{ n*p x k Design matrix (use \code{createX} to make) }
}
\details{
  Let \eqn{\mu_i=X_i beta}, then \eqn{Pr(y_i=j) = exp(\mu_{i,j})/\sum_kexp(\mu_{i,k})}.\cr
  \eqn{X_i} is the submatrix of X corresponding to the
  ith observation.  X has n*p rows.  
  
  Use \code{\link{createX}} to create X. 
}
\value{
  value of log-likelihood (sum of log prob of observed multinomial outcomes).
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{createX}}, \code{\link{rmnlIndepMetrop}} }

\examples{
##
\dontrun{ll=llmnl(beta,y,X)}
}

\keyword{ models }
