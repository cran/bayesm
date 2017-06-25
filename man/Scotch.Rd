\name{Scotch}
\alias{Scotch}
\docType{data}

\title{Survey Data on Brands of Scotch Consumed}

\description{
Data from Simmons Survey. Brands used in last year for those respondents who report consuming scotch. 
}

\usage{data(Scotch)}

\format{
A data frame with 2218 observations on 21 brand variables. \cr
All variables are numeric vectors that are coded 1 if consumed in last year, 0 if not.
}

\source{Edwards, Yancy and Greg Allenby (2003), "Multivariate Analysis of Multiple Response Data," \emph{Journal of Marketing Research}  40, 321--334.}

\references{
Chapter 4, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch\cr \url{http://www.perossi.org/home/bsm-1}}

\examples{
data(Scotch)
cat(" Frequencies of Brands", fill=TRUE)
mat = apply(as.matrix(Scotch), 2, mean)
print(mat)


## use Scotch data to run Multivariate Probit Model
if(0) {
  y = as.matrix(Scotch)
  p = ncol(y)
  n = nrow(y)
  dimnames(y) = NULL
  y = as.vector(t(y))
  y = as.integer(y)
  I_p = diag(p)
  X = rep(I_p,n)
  X = matrix(X, nrow=p)
  X = t(X)
  
  R = 2000
  Data = list(p=p, X=X, y=y)
  Mcmc = list(R=R)
  
  set.seed(66)
  out = rmvpGibbs(Data=Data, Mcmc=Mcmc)
  
  ind = (0:(p-1))*p + (1:p)
  cat(" Betadraws ", fill=TRUE)
  mat = apply(out$betadraw/sqrt(out$sigmadraw[,ind]), 2 , quantile, 
        probs=c(0.01, 0.05, 0.5, 0.95, 0.99))
  attributes(mat)$class = "bayesm.mat"
  summary(mat)
  
  rdraw = matrix(double((R)*p*p), ncol=p*p)
  rdraw = t(apply(out$sigmadraw, 1, nmat))
  attributes(rdraw)$class = "bayesm.var"
  cat(" Draws of Correlation Matrix ", fill=TRUE)
  summary(rdraw)
}
}

\keyword{datasets}
