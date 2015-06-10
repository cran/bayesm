\name{mnlHess}
\alias{mnlHess}
\concept{multinomial logit}
\concept{hessian}


\title{ Computes -Expected Hessian for Multinomial Logit}
\description{
  \code{mnlHess} computes -Expected[Hessian] for Multinomial Logit Model
}
\usage{
mnlHess(beta,y,X)
}
\arguments{
  \item{beta}{ k x 1 vector of coefficients }
  \item{y}{ n x 1 vector of choices, (1, \ldots,p) }
  \item{X}{ n*p x k Design matrix }
}
\details{
  See \code{\link{llmnl}} for information on structure of X array.  Use \code{\link{createX}} to make X.
}
\value{
  k x k matrix
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 3. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{llmnl}}, \code{\link{createX}}, \code{\link{rmnlIndepMetrop}} }
\examples{
##
\dontrun{mnlHess(beta,y,X)}
}
\keyword{ models }
