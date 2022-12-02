\name{llmnl}
\alias{llmnl}
\concept{multinomial logit}
\concept{likelihood}

\title{Evaluate Log Likelihood for Multinomial Logit Model}

\description{\code{llmnl} evaluates log-likelihood for the multinomial logit model.}

\usage{llmnl(beta, y, X)}

\arguments{
  \item{beta}{ \eqn{k x 1} coefficient vector }
  \item{y   }{ \eqn{n x 1} vector of obs on y (1,\ldots, p) }
  \item{X   }{ \eqn{n*p x k} design matrix (use \code{createX} to create \eqn{X}) }
}

\details{
  Let \eqn{\mu_i = X_i beta}, then \eqn{Pr(y_i=j) = exp(\mu_{i,j}) / \sum_k exp(\mu_{i,k})}.\cr
  \eqn{X_i} is the submatrix of \eqn{X} corresponding to the
  \eqn{i}th observation.  \eqn{X} has \eqn{n*p} rows.  
  
  Use \code{\link{createX}} to create \eqn{X}. 
}

\value{Value of log-likelihood (sum of log prob of observed multinomial outcomes).}

\section{Warning}{ This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{createX}}, \code{\link{rmnlIndepMetrop}} }

\examples{
\dontrun{ll=llmnl(beta,y,X)}
}

\keyword{models}
