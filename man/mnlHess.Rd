\name{mnlHess}
\alias{mnlHess}

\title{ Computes -Expected Hessian for Multinomial Logit}
\description{
  \code{mnlHess} computes -Expected[Hessian] for Multinomial Logit Model
}
\usage{
mnlHess(y, X, beta)
}
\arguments{
  \item{y}{ n x 1 vector of choices, (1, \ldots,p) }
  \item{X}{ n*p x k Design matrix }
  \item{beta}{ k x 1 vector of coefficients }
}
\details{
  See \code{\link{llmnl}} for information on structure of X array.  Use \code{\link{createX}} to make X.
}
\value{
  k x k matrix
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 3. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{llmnl}}, \code{\link{createX}}, \code{\link{rmnlIndepMetrop}} }
\examples{
##
\dontrun{mnlHess(y,X,beta)}
}
\keyword{ models }
