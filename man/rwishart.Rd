\name{rwishart}
\alias{rwishart}
\concept{Wishart distribution}
\concept{Inverted Wishart}
\concept{simulation}

\title{ Draw from Wishart and Inverted Wishart Distribution  }

\description{
\code{rwishart} draws from the Wishart and Inverted Wishart distributions.
}

\usage{rwishart(nu, V)}

\arguments{
  \item{nu}{ d.f. parameter}
  \item{V}{ pds location matrix}
}

\details{
  In the parameterization used here, \eqn{W} \eqn{\sim}{~} \eqn{W(nu,V)} with \eqn{E[W]=nuV}.  \cr

  If you want to use an Inverted Wishart prior, you \emph{must invert the location matrix} 
  before calling \code{rwishart}, e.g. \cr
  \eqn{\Sigma} \eqn{\sim}{~} IW(nu ,V);  \eqn{\Sigma^{-1}} \eqn{\sim}{~} \eqn{W(nu, V^{-1})}.
}

\value{
  A list containing:
  \item{W:  }{ Wishart draw }
  \item{IW: }{Inverted Wishart draw}
  \item{C:  }{ Upper tri root of W}
  \item{CI: }{ inv(C), \eqn{W^{-1}} = CICI'}
}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 2, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
set.seed(66)
rwishart(5,diag(3))$IW
}

\keyword{multivariate}
\keyword{distribution}
