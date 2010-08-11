\name{llmnp}
\alias{llmnp}
\concept{multinomial probit}
\concept{GHK method}
\concept{likelihood}

\title{ Evaluate Log Likelihood for Multinomial Probit Model  }
\description{
  \code{llmnp} evaluates the log-likelihood for the multinomial probit model.
}

\usage{
llmnp(beta, Sigma, X, y, r)
}

\arguments{
  \item{beta}{ k x 1 vector of coefficients }
  \item{Sigma}{ (p-1) x (p-1) Covariance matrix of errors }
  \item{X}{ X is n*(p-1) x k array. X is from differenced system. }
  \item{y}{ y is vector of n indicators of multinomial response (1, \ldots, p). }
  \item{r}{ number of draws used in GHK }
}

\details{
  X is (p-1)*n x k matrix.  Use \code{\link{createX}} with \code{DIFF=TRUE} to create X. \cr

  Model for each obs:  \eqn{w = Xbeta + e}.   \eqn{e} \eqn{\sim}{~} \eqn{N(0,Sigma)}. 

   censoring mechanism: 

      if \eqn{y=j (j<p),  w_j > max(w_{-j})} and \eqn{w_j >0}    \cr
      if \eqn{y=p,  w < 0}                   \cr

   To use GHK, we must transform so that these are rectangular regions
   e.g.  if \eqn{y=1,  w_1 > 0} and \eqn{w_1 - w_{-1} > 0}. 

   Define \eqn{A_j} such that if j=1,\ldots,p-1, \eqn{A_jw = A_jmu + A_je > 0} is equivalent to \eqn{y=j}.
   Thus, if y=j, we have \eqn{A_je > -A_jmu}.  Lower truncation is \eqn{-A_jmu} and \eqn{cov = A_jSigmat(A_j)}.  
   For \eqn{j=p}, \eqn{e < - mu}.


}
\value{
  value of log-likelihood (sum of log prob of observed multinomial outcomes).
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapters 2 and 4. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\section{Warning}{
  This routine is a utility routine that does \strong{not} check the
  input arguments for proper dimensions and type.
}

\seealso{ \code{\link{createX}}, \code{\link{rmnpGibbs}} }

\examples{
##
\dontrun{ll=llmnp(beta,Sigma,X,y,r)}
}

\keyword{ models }
