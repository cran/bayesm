\name{rDPGibbs}
\alias{rDPGibbs}
\concept{bayes}
\concept{MCMC}
\concept{normal mixtures}
\concept{Dirichlet Process}
\concept{Gibbs Sampling}

\title{ Density Estimation with Dirichlet Process Prior and Normal Base }

\description{
\code{rDPGibbs} implements a Gibbs Sampler to draw from the posterior for a normal mixture problem with a Dirichlet Process prior. A natural conjugate base prior is used along with priors on the hyper parameters of this distribution. One interpretation of this model is as a normal mixture with a random number of components that can grow with the sample size.
}

\usage{rDPGibbs(Prior, Data, Mcmc)}

\arguments{
  \item{Data }{list(y)}
  \item{Prior}{list(Prioralpha, lambda_hyper)}
  \item{Mcmc }{list(R, keep, nprint, maxuniq, SCALE, gridsize)}
}

\details{
\subsection{Model and Priors}{
  \eqn{y_i} \eqn{\sim}{~} \eqn{N(\mu_i, \Sigma_i)}

  \eqn{\theta_i=(\mu_i,\Sigma_i)} \eqn{\sim}{~} \eqn{DP(G_0(\lambda),alpha)}\cr
  
  \eqn{G_0(\lambda):}\cr
  \eqn{\mu_i | \Sigma_i} \eqn{\sim}{~} \eqn{N(0,\Sigma_i (x) a^{-1})}\cr
  \eqn{\Sigma_i} \eqn{\sim}{~} \eqn{IW(nu,nu*v*I)}
  
  \eqn{\lambda(a,nu,v):}\cr
  \eqn{a} \eqn{\sim}{~} uniform on grid[alim[1], alimb[2]]\cr
  \eqn{nu} \eqn{\sim}{~} uniform on grid[dim(data)-1 + exp(nulim[1]), dim(data)-1 + exp(nulim[2])]\cr
  \eqn{v} \eqn{\sim}{~} uniform on grid[vlim[1], vlim[2]]
  
  \eqn{alpha} \eqn{\sim}{~} \eqn{(1-(\alpha-alphamin)/(alphamax-alphamin))^{power}} \cr
  \eqn{alpha} = alphamin then expected number of components = \code{Istarmin} \cr
  \eqn{alpha} = alphamax then expected number of components = \code{Istarmax} 

  We parameterize the prior on \eqn{\Sigma_i} such that \eqn{mode(\Sigma)= nu/(nu+2) vI}. The support of nu enforces valid IW density; \eqn{nulim[1] > 0}

  We use the structure for \code{nmix} that is compatible with the \code{bayesm} routines for finite mixtures of normals. This allows us to use the same summary and plotting methods.  

  The default choices of \code{alim}, \code{nulim}, and \code{vlim} determine the location and approximate size of candidate "atoms" or possible normal components. The defaults are sensible given that we scale the data.  Without scaling, you want to insure that \code{alim} is set for a wide enough range of values (remember a is a precision parameter) and the \code{v} is big enough to propose \code{Sigma} matrices wide enough to cover the data range.  

  A careful analyst should look at the posterior distribution of \code{a}, \code{nu}, \code{v} to make sure that the support is set correctly in \code{alim}, \code{nulim}, \code{vlim}.  In other words, if we see the posterior bunched up at one end of these support ranges, we should widen the range and rerun.  

  If you want to force the procedure to use many small atoms, then set \code{nulim} to consider only large values and set \code{vlim} to consider only small scaling constants. Set \code{Istarmax} to a large number.  This will create a very "lumpy" density estimate somewhat like the classical Kernel density estimates. Of course, this is not advised if you have a prior belief that densities are relatively smooth.
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(y)}}
  \tabular{ll}{
    \code{y: } \tab \eqn{n x k} matrix of observations on \eqn{k} dimensional data
    }
  \emph{\code{Prior = list(Prioralpha, lambda_hyper)} [optional]}
  \tabular{ll}{
    \code{Prioralpha:  } \tab \code{list(Istarmin, Istarmax, power)} \cr
      \code{$Istarmin:  } \tab is expected number of components at lower bound of support of alpha (def: 1) \cr
      \code{$Istarmax:  } \tab is expected number of components at upper bound of support of alpha (def: \code{min(50, 0.1*nrow(y))}) \cr
      \code{$power:     } \tab is the power parameter for alpha prior (def: 0.8) \cr
    \code{lambda_hyper:} \tab \code{list(alim, nulim, vlim)} \cr
      \code{$alim:      } \tab defines support of a distribution (def: \code{c(0.01, 10)}) \cr
      \code{$nulim:     } \tab defines support of nu distribution (def: \code{c(0.01, 3)}) \cr
      \code{$vlim:      } \tab defines support of v distribution (def: \code{c(0.1, 4)})
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, maxuniq, SCALE, gridsize)} [only \code{R} required]}
  \tabular{ll}{
   \code{R:        } \tab number of MCMC draws \cr
   \code{keep:     } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
   \code{nprint:   } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
   \code{maxuniq:  } \tab storage constraint on the number of unique components (def: 200) \cr
   \code{SCALE:    } \tab should data be scaled by mean,std deviation before posterior draws (def: \code{TRUE}) \cr
   \code{gridsize: } \tab number of discrete points for hyperparameter priors (def: 20)
    }
}
\subsection{\code{nmix} Details}{
  \code{nmix} is a list with 3 components. Several functions in the \code{bayesm} package that involve a Dirichlet Process or mixture-of-normals return \code{nmix}. Across these functions, a common structure is used for \code{nmix} in order to utilize generic summary and plotting functions. 
  \tabular{ll}{
  \code{probdraw:} \tab \eqn{ncomp x R/keep} matrix that reports the probability that each draw came from a particular component \cr
  \code{zdraw:   } \tab \eqn{R/keep x nobs} matrix that indicates which component each draw is assigned to \cr
  \code{compdraw:} \tab A list of \eqn{R/keep} lists of \eqn{ncomp} lists. Each of the inner-most lists has 2 elemens: a vector of draws for \code{mu} and a matrix of draws for the Cholesky root of \code{Sigma}.
  }
}
}

\value{
A list containing:
  \item{nmix      }{ a list containing: \code{probdraw}, \code{zdraw}, \code{compdraw} (see \dQuote{\code{nmix} Details} section)}
  \item{alphadraw }{ \eqn{R/keep x 1} vector of alpha draws}
  \item{nudraw    }{ \eqn{R/keep x 1} vector of nu draws}
  \item{adraw     }{ \eqn{R/keep x 1} vector of a draws}
  \item{vdraw     }{ \eqn{R/keep x 1} vector of v draws}
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\seealso{ \code{\link{rnmixGibbs}}, \code{\link{rmixture}}, \code{\link{rmixGibbs}}, 
          \code{\link{eMixMargDen}}, \code{\link{momMix}}, \code{\link{mixDen}}, \code{\link{mixDenBi}}}

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)

## simulate univariate data from Chi-Sq

N = 200
chisqdf = 8
y1 = as.matrix(rchisq(N,df=chisqdf))

## set arguments for rDPGibbs

Data1 = list(y=y1)
Prioralpha = list(Istarmin=1, Istarmax=10, power=0.8)
Prior1 = list(Prioralpha=Prioralpha)
Mcmc = list(R=R, keep=1, maxuniq=200)

out1 = rDPGibbs(Prior=Prior1, Data=Data1, Mcmc=Mcmc)

if(0){
  ## plotting examples
  rgi = c(0,20)
  grid = matrix(seq(from=rgi[1],to=rgi[2],length.out=50), ncol=1)
  deltax = (rgi[2]-rgi[1]) / nrow(grid)
  plot(out1$nmix, Grid=grid, Data=y1)
  
  ## plot true density with historgram
  plot(range(grid[,1]), 1.5*range(dchisq(grid[,1],df=chisqdf)),
       type="n", xlab=paste("Chisq ; ",N," obs",sep=""), ylab="")
  hist(y1, xlim=rgi, freq=FALSE, col="yellow", breaks=20, add=TRUE)
  lines(grid[,1], dchisq(grid[,1],df=chisqdf) / (sum(dchisq(grid[,1],df=chisqdf))*deltax),
        col="blue", lwd=2)
}

## simulate bivariate data from the  "Banana" distribution (Meng and Barnard) 

banana = function(A, B, C1, C2, N, keep=10, init=10) { 
  R = init*keep + N*keep
  x1 = x2 = 0
  bimat = matrix(double(2*N), ncol=2)
  for (r in 1:R) { 
    x1 = rnorm(1,mean=(B*x2+C1) / (A*(x2^2)+1), sd=sqrt(1/(A*(x2^2)+1)))
    x2 = rnorm(1,mean=(B*x2+C2) / (A*(x1^2)+1), sd=sqrt(1/(A*(x1^2)+1)))
    if (r>init*keep && r\%\%keep==0) {
      mkeep = r/keep
      bimat[mkeep-init,] = c(x1,x2)
    } 
  }
  return(bimat)
}

set.seed(66)
nvar2 = 2
A = 0.5
B = 0
C1 = C2 = 3
y2 = banana(A=A, B=B, C1=C1, C2=C2, 1000)

Data2 = list(y=y2)
Prioralpha = list(Istarmin=1, Istarmax=10, power=0.8)
Prior2 = list(Prioralpha=Prioralpha)
Mcmc = list(R=R, keep=1, maxuniq=200)

out2 = rDPGibbs(Prior=Prior2, Data=Data2, Mcmc=Mcmc)


if(0){
  ## plotting examples
  
  rx1 = range(y2[,1])
  rx2 = range(y2[,2])
  x1 = seq(from=rx1[1], to=rx1[2], length.out=50)
  x2 = seq(from=rx2[1], to=rx2[2], length.out=50)
  grid = cbind(x1,x2)
  plot(out2$nmix, Grid=grid, Data=y2)
  
  ## plot true bivariate density
  tden = matrix(double(50*50), ncol=50)
  for (i in 1:50) { 
  for (j in 1:50) {
        tden[i,j] = exp(-0.5*(A*(x1[i]^2)*(x2[j]^2) + 
                    (x1[i]^2) + (x2[j]^2) - 2*B*x1[i]*x2[j] - 
                    2*C1*x1[i] - 2*C2*x2[j]))
  }}
  tden = tden / sum(tden)
  image(x1, x2, tden, col=terrain.colors(100), xlab="", ylab="")
  contour(x1, x2, tden, add=TRUE, drawlabels=FALSE)
  title("True Density")
}
}

\keyword{multivariate}
