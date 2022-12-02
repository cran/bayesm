\name{breg}
\alias{breg}
\concept{bayes}
\concept{regression}

\title{Posterior Draws from a Univariate Regression with Unit Error Variance}

\description{
  \code{breg} makes one draw from the posterior of a univariate regression
  (scalar dependent variable) given the error variance = 1.0.
  A natural conjugate (normal) prior is used.
}

\usage{breg(y, X, betabar, A)}

\arguments{
  \item{y       }{ \eqn{n x 1} vector of values of dep variable }
  \item{X       }{ \eqn{n x k} design matrix}
  \item{betabar }{ \eqn{k x 1} vector for the prior mean of the regression coefficients }
  \item{A       }{ \eqn{k x k} prior precision matrix }
}

\details{
  model: \eqn{y = X'\beta + e} with   \eqn{e} \eqn{\sim}{~} \eqn{N(0,1)}.  \cr
  prior: \eqn{\beta} \eqn{\sim}{~} \eqn{N(betabar, A^{-1})}.
}

\value{\eqn{k x 1} vector containing a draw from the posterior distribution}

\section{Warning}{
This routine is a utility routine that does \strong{not} check theinput arguments for proper dimensions and type. In particular, \code{X} must be a matrix. If you have a vector for \code{X}, coerce itinto a matrix with one column.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\references{For further discussion, see \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch. }

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=1000} else {R=10}

## simulate data
set.seed(66)
n = 100
X = cbind(rep(1,n), runif(n)); beta = c(1,2)
y = X \%*\% beta + rnorm(n)

## set prior
betabar = c(0,0)
A = diag(c(0.05, 0.05))

## make draws from posterior
betadraw = matrix(double(R*2), ncol = 2)
for (rep in 1:R) {betadraw[rep,] = breg(y,X,betabar,A)}

## summarize draws
mat = apply(betadraw, 2, quantile, probs=c(0.01, 0.05, 0.50, 0.95, 0.99))
mat = rbind(beta,mat); rownames(mat)[1] = "beta"
print(mat)
}

\keyword{models}
\keyword{regression}
\keyword{distribution}
