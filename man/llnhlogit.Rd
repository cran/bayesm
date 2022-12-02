\name{llnhlogit}
\alias{llnhlogit}
\concept{multinomial logit}
\concept{non-homothetic utility}

\title{Evaluate Log Likelihood for non-homothetic Logit Model}

\description{
\code{llnhlogit} evaluates log-likelihood for the Non-homothetic Logit model.
}

\usage{llnhlogit(theta, choice, lnprices, Xexpend)}

\arguments{
  \item{theta   }{ parameter vector (see details section) }
  \item{choice  }{ \eqn{n x 1} vector of choice (1,\ldots,p) }
  \item{lnprices}{ \eqn{n x p}  array of log-prices}
  \item{Xexpend }{ \eqn{n x d} array of vars predicting expenditure }
}

\details{
  Non-homothetic logit model, \eqn{Pr(i) = exp(tau v_i) / sum_j exp(tau v_j)} \cr
  
  \eqn{v_i = alpha_i - e^{kappaStar_i}u^i - lnp_i} \cr
  tau is the scale parameter of extreme value error distribution.\cr
  \eqn{u^i} solves \eqn{u^i = psi_i(u^i)E/p_i}.\cr
  \eqn{ln(psi_i(U)) = alpha_i - e^{kappaStar_i}U}. \cr
  \eqn{ln(E) = gamma'Xexpend}.\cr

  Structure of theta vector: \cr
  alpha: \eqn{p x 1} vector of utility intercepts.\cr
  kappaStar: \eqn{p x 1} vector of utility rotation parms expressed on natural log scale. \cr
  gamma: \eqn{k x 1} -- expenditure variable coefs.\cr
  tau: \eqn{1 x 1} -- logit scale parameter.\cr
}

\value{Value of log-likelihood (sum of log prob of observed multinomial outcomes).}

\section{Warning}{
This routine is a utility routine that does \strong{not} check the input arguments for proper dimensions and type.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see Chapter 4, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\seealso{\code{\link{simnhlogit}}}

\examples{
N=1000; p=3; k=1
theta = c(rep(1,p), seq(from=-1,to=1,length=p), rep(2,k), 0.5)
lnprices = matrix(runif(N*p), ncol=p)
Xexpend = matrix(runif(N*k), ncol=k)
simdata = simnhlogit(theta, lnprices, Xexpend)

## evaluate likelihood at true theta
llstar = llnhlogit(theta, simdata$y, simdata$lnprices, simdata$Xexpend)
}

\keyword{models}