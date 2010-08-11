\name{llnhlogit}
\alias{llnhlogit}
\concept{multinomial logit}
\concept{non-homothetic utility}

\title{ Evaluate Log Likelihood for non-homothetic Logit Model  }
\description{
  \code{llmnp} evaluates log-likelihood for the Non-homothetic Logit model.
}

\usage{
llnhlogit(theta, choice, lnprices, Xexpend)
}

\arguments{
  \item{theta}{ parameter vector (see details section) }
  \item{choice}{ n x 1 vector of choice (1, \ldots, p) }
  \item{lnprices}{ n x p  array of log-prices}
  \item{Xexpend}{ n x d array of vars predicting expenditure }
}

\details{
  Non-homothetic logit model with: \eqn{ln(psi_i(U)) = alpha_i - e^{k_i}U} \cr

  Structure of theta vector \cr
  alpha: (p x 1) vector of utility intercepts.\cr
  k: (p x 1) vector of utility rotation parms. \cr
  gamma: (k x 1) -- expenditure variable coefs.\cr
  tau: (1 x 1) -- logit scale parameter.\cr
}

\value{
  value of log-likelihood (sum of log prob of observed multinomial outcomes).
}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch,Chapter 4. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{simnhlogit}} }
\examples{
##
\dontrun{ll=llnhlogit(theta,choice,lnprices,Xexpend)}
}
\keyword{ models }
