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
  \item{Prior}{ list(nu,V,mubar,Am,gsigma,gl11,gl22,gl12,Lambdanu,LambdaV,ge) optional }
  \item{Mcmc}{ list(R,keep,ndghk,printevery,e,y,mu,Sigma,sigma,tau,Lambda) optional }
}
\details{
  Model: n=nrow(X) individuals respond to m=ncol(X) questions. all questions are on a scale 1, \ldots, k.\cr

  for respondent i and question j, \cr
  \eqn{x_ij = d}, if \eqn{c_d-1 \le y_i,j \le c_d}. \cr
  d=1,\ldots,k. \eqn{c_d = a + bd +ed^2}. \cr

  \eqn{y_i = mu + tau_i*iota + sigma_i*z_i}.  \eqn{z_i} \eqn{\sim}{~} \eqn{N(0,Sigma)}. \cr

  \eqn{(tau_i,ln(sigma_i))} \eqn{\sim}{~} \eqn{N(phi,Lamda)}.  \eqn{phi=(0,lambda_{22})}. \cr

  Priors:\cr
  mu \eqn{\sim}{~} \eqn{N(mubar, Am{^-1})}.\cr
  Sigma \eqn{\sim}{~} IW(nu,V).\cr
  Lambda \eqn{\sim}{~} IW(Lambdanu,LambdaV).\cr
  e \eqn{\sim}{~} unif on a grid. \cr
}
\value{
  a list containing:
  \item{drawSigma}{R/keep x m*m array of Sigma draws}
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

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Allenby, McCulloch, and Rossi, Case Study on Scale Usage Heterogeneity. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Rob McCulloch and Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}
\keyword{ models }
