\name{rwishart}
\alias{rwishart}
\title{ Draw from Wishart and Inverted Wishart Distribution  }
\description{
  \code{rwishart} draws from the Wishart and Inverted Wishart distributions.
}
\usage{
rwishart(nu, V)
}
\arguments{
  \item{nu}{ d.f. parameter}
  \item{V}{ pds location matrix}
}
\details{
  In the parameterization used here, eqn{W} \eqn{\sim}{~} \eqn{W(nu,V)}, \eqn{E[W]=nuV}.  \cr

  If you want to use an Inverted Wishart prior, you \emph{must invert the location matrix} 
  before calling \code{rwishart}, e.g. \cr
  \eqn{Sigma} \eqn{\sim}{~} IW(nu,V);  \eqn{Sigma^{-1}} \eqn{\sim}{~} \eqn{W(nu,V^{-1})}.
}
\value{
  \item{W}{ Wishart draw }
  \item{IW }{Inverted Wishart draw}
  \item{C }{ Upper tri root of W}
  \item{CI }{ inv(C), \eqn{W^{-1}} = CICI'}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Chapter 2. \cr
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
rwishart(5,diag(3))$IW
}
\keyword{ multivariate }
