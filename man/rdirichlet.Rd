\name{rdirichlet}
\alias{rdirichlet}
\concept{dirichlet distribution}
\concept{simulation}

\title{Draw From Dirichlet Distribution}

\description{
\code{rdirichlet} draws from Dirichlet
}

\usage{rdirichlet(alpha)}

\arguments{
  \item{alpha}{vector of Dirichlet parms (must be > 0)}
}

\value{Vector of draws from Dirichlet}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 2, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
set.seed(66)
rdirichlet(c(rep(3,5)))
}

\keyword{distribution}

