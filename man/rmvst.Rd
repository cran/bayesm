\name{rmvst}
\alias{rmvst}
\concept{multivariate t distribution}
\concept{student-t}
\concept{simulation}

\title{ Draw from Multivariate Student-t }

\description{
  \code{rmvst} draws from a multivariate student-t distribution.
}

\usage{rmvst(nu, mu, root)}

\arguments{
  \item{nu}{ d.f. parameter }
  \item{mu}{ mean vector }
  \item{root}{ Upper Tri Cholesky Root of Sigma }
}

\value{length(mu) draw vector}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see \emph{Bayesian Statistics and Marketing}  by Rossi, Allenby, and McCulloch.}

\seealso{ \code{\link{lndMvst}}}

\examples{
set.seed(66)
rmvst(nu=5, mu=c(rep(0,2)), root=chol(matrix(c(2,1,1,2), ncol=2)))
}

\keyword{multivariate}
\keyword{distribution}
