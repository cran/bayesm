\name{createX}
\alias{createX}
\concept{multinomial logit}
\concept{multinomial probit}

\title{Create X Matrix for Use in Multinomial Logit and Probit Routines}

\description{
  \code{createX} makes up an X matrix in the form expected by Multinomial
  Logit (\code{\link{rmnlIndepMetrop}} and \code{\link{rhierMnlRwMixture}})
  and Probit (\code{\link{rmnpGibbs}} and \code{\link{rmvpGibbs}}) routines.  
  Requires an array of alternative-specific variables and/or an
  array of "demographics" (or variables constant across alternatives) which
  may vary across choice occasions. 
}
\usage{createX(p, na, nd, Xa, Xd, INT = TRUE, DIFF = FALSE, base=p)}

\arguments{
  \item{p}{ integer number of choice alternatives }
  \item{na}{ integer number of alternative-specific vars in \code{Xa}  }
  \item{nd}{ integer number of non-alternative specific vars }
  \item{Xa}{ \eqn{n x p*na} matrix of alternative-specific vars }
  \item{Xd}{ \eqn{n x nd} matrix of non-alternative specific vars }
  \item{INT}{ logical flag for inclusion of intercepts }
  \item{DIFF}{ logical flag for differencing wrt to base alternative }
  \item{base}{ integer index of base choice alternative }
  Note: \code{na}, \code{nd}, \code{Xa}, \code{Xd} can be \code{NULL} to indicate lack of \code{Xa} or \code{Xd} variables.
}

\value{\code{X} matrix of dimension \eqn{n*(p-DIFF) x [(INT+nd)*(p-1) + na]}.}

\note{\code{\link{rmnpGibbs}} assumes that the \code{base} alternative is the default.}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\seealso{\code{\link{rmnlIndepMetrop}}, \code{\link{rmnpGibbs}} }

\examples{
na=2; nd=1; p=3
vec = c(1, 1.5, 0.5, 2, 3, 1, 3, 4.5, 1.5)
Xa = matrix(vec, byrow=TRUE, ncol=3)
Xa = cbind(Xa,-Xa)
Xd = matrix(c(-1,-2,-3), ncol=1)
createX(p=p, na=na, nd=nd, Xa=Xa, Xd=Xd)
createX(p=p, na=na, nd=nd, Xa=Xa, Xd=Xd, base=1)
createX(p=p, na=na, nd=nd, Xa=Xa, Xd=Xd, DIFF=TRUE)
createX(p=p, na=na, nd=nd, Xa=Xa, Xd=Xd, DIFF=TRUE, base=2)
createX(p=p, na=na, nd=NULL, Xa=Xa, Xd=NULL)
createX(p=p, na=NULL, nd=nd, Xa=NULL, Xd=Xd)
}

\keyword{array}
\keyword{utilities}
