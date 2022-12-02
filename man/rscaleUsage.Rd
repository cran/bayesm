\name{rscaleUsage}
\alias{rscaleUsage}
\concept{MCMC}
\concept{bayes}
\concept{ordinal data}
\concept{scale usage}
\concept{hierarchical models}

\title{MCMC Algorithm for Multivariate Ordinal Data with Scale Usage Heterogeneity}

\description{
\code{rscaleUsage} implements an MCMC algorithm for multivariate ordinal data with scale usage heterogeniety.  
}

\usage{rscaleUsage(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(x, k)}
  \item{Prior}{list(nu, V, mubar, Am, gs, Lambdanu, LambdaV)}
  \item{Mcmc }{list(R, keep, nprint, ndghk, e, y, mu, Sigma, sigma, tau, Lambda)}
}

\details{
\subsection{Model and Priors}{
  \eqn{n} = \code{nrow(x)} individuals respond to \eqn{p} = \code{ncol(x)} questions; 
  all questions are on a scale \eqn{1, \ldots, k} for respondent \eqn{i} and question \eqn{j}, \cr
  
  \eqn{x_{ij} = d} if \eqn{c_{d-1} \le y_{ij} \le c_d} where \eqn{d = 1, \ldots, k} and \eqn{c_d = a + bd + ed^2} \cr

  \eqn{y_i = mu + tau_i*iota + sigma_i*z_i} with \eqn{z_i} \eqn{\sim}{~} \eqn{N(0, Sigma)}

  \eqn{(tau_i, ln(sigma_i))} \eqn{\sim}{~} \eqn{N(\phi, Lamda)}\cr
  \eqn{\phi = (0, lambda_{22})} \cr
  \eqn{mu} \eqn{\sim}{~} \eqn{N(mubar, Am^{-1})}\cr
  \eqn{Sigma} \eqn{\sim}{~} \eqn{IW(nu, V)} \cr
  \eqn{Lambda} \eqn{\sim}{~} \eqn{IW(Lambdanu, LambdaV)} \cr
  \eqn{e} \eqn{\sim}{~} unif on a grid 
  
  It is highly recommended that the user choose the default prior settings. 
  If you wish to change prior settings and/or the grids used, please carefully read the case study listed in the reference below.
}
\subsection{Argument Details}{
  \emph{\code{Data  = list(x, k)}}
  \tabular{ll}{
    \code{x: } \tab \eqn{n x p} matrix of discrete responses \cr
    \code{k: } \tab number of discrete rating scale options
    }
  \emph{\code{Prior = list(nu, V, mubar, Am, gs, Lambdanu, LambdaV)} [optional]}
  \tabular{ll}{
    \code{nu:       } \tab d.f. parameter for Sigma prior (def: p + 3) \cr
    \code{V:        } \tab scale location matrix for Sigma prior (def: nu*I) \cr
    \code{mubar:    } \tab \eqn{p x 1} vector of prior means (def: \code{rep(k/2,p)}) \cr
    \code{Am:       } \tab \eqn{p x p} prior precision matrix (def: 0.01*I) \cr
    \code{gs:       } \tab grid size for sigma (def: 100) \cr
    \code{Lambdanu: } \tab d.f. parameter for Lambda prior (def: 20) \cr
    \code{LambdaV:  } \tab scale location matrix for Lambda prior (def: (Lambdanu - 3)*Lambda)
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, ndghk, e, y, mu, Sigma, sigma, tau, Lambda)} [only \code{R} required]}
  \tabular{ll}{
    \code{R:      } \tab number of MCMC draws (def: 1000) \cr
    \code{keep:   } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint: } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{ndghk:  } \tab number of draws for a GHK integration (def: 100) \cr
    \code{e:      } \tab initial value (def: 0) \cr
    \code{y:      } \tab initial values (def: x) \cr
    \code{mu:     } \tab initial values (def: \code{apply(y,2,mean)}, a p-length vector) \cr
    \code{Sigma:  } \tab initial value (def: \code{var(y)}) \cr
    \code{sigma:  } \tab initial values (def: \code{rep(1,n)}) \cr
    \code{tau:    } \tab initial values (def: \code{rep(0,n)}) \cr
    \code{Lambda: } \tab initial values (def: \code{matrix(c(4,0,0,.5),ncol=2)})
    }
}
}

\value{
  A list containing:
  \item{Sigmadraw  }{\eqn{R/keep x p*p} matrix of Sigma draws -- each row is the vector form of Sigma}
  \item{mudraw     }{\eqn{R/keep x p} matrix of mu draws}
  \item{taudraw    }{\eqn{R/keep x n} matrix of tau draws}
  \item{sigmadraw  }{\eqn{R/keep x n} matrix of sigma draws}
  \item{Lambdadraw }{\eqn{R/keep x 4} matrix of Lamda draws}
  \item{edraw      }{\eqn{R/keep x 1} vector of e draws}
}

\section{Warning}{
\eqn{tau_i}, \eqn{sigma_i} are identified from the scale usage patterns in the \eqn{p} questions asked per respondent (# cols of \eqn{x}).  Do not attempt to use this on datasets with only a small number of total questions.
}

\author{Rob McCulloch (Arizona State University) and Peter Rossi (Anderson School, UCLA), \email{perossichi@gmail.com}.}

\references{ For further discussion, see Case Study 3 on Overcoming Scale Usage Heterogeneity, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=1000} else {R=5} 
set.seed(66)

data(customerSat)
surveydat = list(k=10, x=as.matrix(customerSat))

out = rscaleUsage(Data=surveydat, Mcmc=list(R=R))
summary(out$mudraw)
}

\keyword{models}
