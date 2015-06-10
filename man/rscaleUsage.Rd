\name{rscaleUsage}
\alias{rscaleUsage}
\concept{MCMC}
\concept{bayes}
\concept{ordinal data}
\concept{scale usage}
\concept{hierarchical models}

\title{ MCMC Algorithm for Multivariate Ordinal Data with Scale Usage Heterogeneity.}
\description{
  \code{rscaleUsage} implements an MCMC algorithm for multivariate ordinal data with scale usage heterogeniety.  
}
\usage{
rscaleUsage(Data,Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(k,x)}
  \item{Prior}{ list(nu,V,mubar,Am,gsigma,gl11,gl22,gl12,Lambdanu,LambdaV,ge) (optional) }
  \item{Mcmc}{ list(R,keep,ndghk,nprint,e,y,mu,Sigma,sigma,tau,Lambda) (optional) }
}
\details{
  Model: n=nrow(x) individuals respond to m=ncol(x) questions. all questions are on a scale 1, \ldots, k.
  for respondent i and question j, \cr
  \eqn{x_{ij} = d}, if \eqn{c_{d-1} \le y_{ij} \le c_d}. \cr
  d=1,\ldots,k. \eqn{c_d = a + bd +ed^2}. \cr

  \eqn{y_i = mu + tau_i*iota + sigma_i*z_i}.  \eqn{z_i} \eqn{\sim}{~} \eqn{N(0,Sigma)}. \cr



  Priors:\cr
  \eqn{(tau_i,ln(sigma_i))} \eqn{\sim}{~} \eqn{N(\phi,Lamda)}.  \eqn{\phi=(0,lambda_{22})}. \cr
  mu \eqn{\sim}{~} \eqn{N(mubar, Am{^-1})}.\cr
  Sigma \eqn{\sim}{~} IW(nu,V).\cr
  Lambda \eqn{\sim}{~} IW(Lambdanu,LambdaV).\cr
  e \eqn{\sim}{~} unif on a grid. \cr
}
\value{
  a list containing:
  \item{Sigmadraw}{R/keep x m*m array of Sigma draws}
  \item{mudraw}{R/keep x m array of mu draws}
  \item{taudraw}{R/keep x n array of tau draws}
  \item{sigmadraw}{R/keep x n array of sigma draws}
  \item{Lambdadraw}{R/keep x 4 array of Lamda draws}
  \item{edraw}{R/keep x  1 array of  e draws}
}
\note{
  It is \strong{highly} recommended that the user choose the default settings. This means not specifying the argument
  \code{Prior} and setting \code{R} in Mcmc and \code{Data} only.  If you wish to change prior settings and/or
  the grids used, please read the case study in Allenby et al carefully.
}
\section{Warning}{
  \eqn{tau_i}, \eqn{sigma_i} are identified from the scale usage patterns in the m questions asked per
   respondent (\# cols of x).  Do not attempt to use this on data sets with only a small number of
   total questions!
  
}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby, and McCulloch, Case Study on Scale Usage Heterogeneity. \cr
  \url{http://www.perossi.org/home/bsm-1}
}

\author{ Rob McCulloch and Peter Rossi, Anderson School, UCLA,
  \email{perossichi@gmail.com}.
}

\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=1000} else {R=1} 
{
data(customerSat)
surveydat = list(k=10,x=as.matrix(customerSat))

Mcmc1 = list(R=R)
set.seed(66)
out=rscaleUsage(Data=surveydat,Mcmc=Mcmc1)

summary(out$mudraw)

}
}
\keyword{ models }
