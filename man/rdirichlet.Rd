\name{rdirichlet}
\alias{rdirichlet}
\concept{dirichlet distribution}
\concept{simulation}

\title{ Draw From Dirichlet Distribution }
\description{
  \code{rdirichlet} draws from Dirichlet
}
\usage{
rdirichlet(alpha)
}
\arguments{
  \item{alpha}{ vector of Dirichlet parms (must be > 0)}
}

\value{
  Vector of draws from Dirichlet
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

\examples{
##
set.seed(66)
rdirichlet(c(rep(3,5)))
}
\keyword{ distribution }

