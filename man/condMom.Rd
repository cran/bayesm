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
  \item{sigi}{ length(x) dim inverse of covariance matrix }
  \item{i}{ conditional distribution of ith element }
}
\details{
  \eqn{x} \eqn{\sim}{~} \eqn{MVN(mu,sigi^{-1})}.

  \code{condMom} computes moments of \eqn{x_i} given \eqn{x_{-i}}.
}
\value{
  a list containing:

  \item{cmean }{ cond mean }
  \item{cvar }{ cond variance}
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
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
