\name{mnlHess}
\alias{mnlHess}
\concept{multinomial logit}
\concept{hessian}


\title{ Computes --Expected Hessian for Multinomial Logit}

\description{\code{mnlHess} computes expected Hessian (\eqn{E[H]}) for Multinomial Logit Model.}

\usage{mnlHess(beta, y, X)}

\arguments{
  \item{beta}{ \eqn{k x 1} vector of coefficients }
  \item{y}{ \eqn{n x 1} vector of choices, (\eqn{1,\ldots,p}) }
  \item{X}{ \eqn{n*p x k} Design matrix }
}
\details{
  See \code{\link{llmnl}} for information on structure of X array.  Use \code{\link{createX}} to make X.
}

\value{\eqn{k x k} matrix}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\seealso{ \code{\link{llmnl}}, \code{\link{createX}}, \code{\link{rmnlIndepMetrop}} }

\examples{
\dontrun{mnlHess(beta, y, X)}
}

\keyword{models}
