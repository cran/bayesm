\name{condMom}
\alias{condMom}
\concept{normal distribution}
\concept{conditional distribution}

\title{ Computes Conditional Mean/Var of One Element of MVN given All Others }
\description{
  \code{condMom} compute moments of conditional distribution of ith element of normal given
 all others.
}
\usage{
condMom(x, mu, sigi, i)
}
\arguments{
  \item{x}{ vector of values to condition on - ith element not used }
  \item{mu}{ length(x) mean vector }
  \item{sigi}{ length(x)-dim covariance matrix }
  \item{i}{ conditional distribution of ith element }
}
\details{
  \eqn{x} \eqn{\sim}{~} \eqn{MVN(mu,Sigma)}.

  \code{condMom} computes moments of \eqn{x_i} given \eqn{x_{-i}}.
}
\value{
  a list containing:

  \item{cmean }{ cond mean }
  \item{cvar }{ cond variance}
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
sig=matrix(c(1,.5,.5,.5,1,.5,.5,.5,1),ncol=3)
sigi=chol2inv(chol(sig))
mu=c(1,2,3)
x=c(1,1,1)
condMom(x,mu,sigi,2)

}
\keyword{ distribution }
