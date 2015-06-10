\name{ghkvec}
\alias{ghkvec}
\concept{multivariate normal distribution}
\concept{GHK method}
\concept{integral}

\title{ Compute GHK approximation to Multivariate Normal Integrals }
\description{
  \code{ghkvec} computes the GHK approximation to the integral of a
  multivariate normal density over a half plane defined by a set
  of truncation points.  
}
\usage{
ghkvec(L, trunpt, above, r, HALTON=TRUE, pn)
}
\arguments{
  \item{L}{ lower triangular Cholesky root of covariance matrix }
  \item{trunpt}{ vector of truncation points}
  \item{above}{ vector of indicators for truncation above(1) or below(0) }
  \item{r}{ number of draws to use in GHK }
  \item{HALTON}{ if TRUE, use Halton sequence. If FALSE, use R::runif random number generator (optional / def: TRUE)}
  \item{pn}{ prime number used for Halton sequence (optional / def: the smallest prime numbers, i.e. 2, 3, 5, ...)}  
}
\value{
  approximation to integral
}
\note{
  \code{ghkvec} can accept a vector of truncations and compute more than one
  integral.  That is, length(trunpt)/length(above) number of different integrals,
  each with the same Sigma and mean 0 but different truncation points. See 
  example below for an example with two integrals at different truncation points. \cr
  
  User can choose what random number to use for the numerical integration: psuedo-random numbers by \code{R::runif} or quasi-random numbers by Halton sequence. Generally, the quasi-random sequence (e.g., Halton) is more uniformly distributed within domain, so it shows lower error and improved convergence than the psuedo-random sequence (Morokoff and Caflisch, 1995). \cr 
  
  For the prime numbers generating Halton sequence, we suggest to use the first smallest prime numbers. Halton (1960) and Kocis and Whiten (1997) prove that their discrepancy measures (how uniformly the sample points are distributed) have the upper bounds, which decrease as the generating prime number decreases. \cr
  
  Note: For a high dimensional integration (10 or more dimension), we suggest to use the psuedo-random number generator (\code{R::runif}). According to Kocis and Whiten (1997), Halton sequences may be highly correlated when the dimension is 10 or more.
}
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch,  Chapter 2. \cr
  \url{http://www.perossi.org/home/bsm-1} \cr
  
  For Halton sequence, see Halton (1960, Numerische Mathematik), Morokoff and Caflisch (1995, Journal of Computational Physics), and Kocis and Whiten (1997, ACM Transactions on Mathematical Software).
}

\author{
  Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.\cr
  Keunwoo Kim, Anderson School, UCLA, \email{keunwoo.kim@gmail.com}
}

\examples{
Sigma=matrix(c(1,.5,.5,1),ncol=2)
L=t(chol(Sigma))
trunpt=c(0,0,1,1)
above=c(1,1)

# drawn by Halton sequence
ghkvec(L,trunpt,above,r=100)

# use prime number 11 and 13
ghkvec(L,trunpt,above,r=100,HALTON=TRUE,pn=c(11,13))

# drawn by R::runif
ghkvec(L,trunpt,above,r=100,HALTON=FALSE)
}
\keyword{ distribution }
