\name{rbayesBLP}
\alias{rbayesBLP}
\concept{bayes}
\concept{random coefficient logit}
\concept{BLP}
\concept{Metropolis Hasting}

\title{Bayesian Analysis of Random Coefficient Logit Models Using Aggregate Data}

\description{
\code{rbayesBLP} implements a hybrid MCMC algorithm for aggregate level sales data in a market with differentiated products. bayesm version 3.1-0 and prior verions contain an error when using instruments with this function; this will be fixed in a future version.
}

\usage{rbayesBLP(Data, Prior, Mcmc)}

\arguments{
  \item{Data }{list(X, share, J, Z)}
  \item{Prior}{list(sigmasqR, theta_hat, A, deltabar, Ad, nu0, s0_sq, VOmega)}
  \item{Mcmc }{list(R, keep, nprint, H, initial_theta_bar, initial_r, initial_tau_sq, initial_Omega, initial_delta, s, cand_cov, tol)}
}

\value{
  A list containing:
  \item{thetabardraw }{\eqn{K x R/keep} matrix of random coefficient mean draws}
  \item{Sigmadraw    }{\eqn{K*K x R/keep} matrix of random coefficient variance draws}
  \item{rdraw        }{\eqn{K*K x R/keep} matrix of \eqn{r} draws (same information as in \code{Sigmadraw})}
  \item{tausqdraw    }{\eqn{R/keep x 1} vector of aggregate demand shock variance draws}
  \item{Omegadraw    }{\eqn{2*2 x R/keep} matrix of correlated endogenous shock variance draws}
  \item{deltadraw    }{\eqn{I x R/keep} matrix of endogenous structural equation coefficient draws}
  \item{acceptrate   }{scalor of acceptance rate of Metropolis-Hasting}
  \item{s            }{scale parameter used for Metropolis-Hasting}
  \item{cand_cov     }{var-cov matrix used for Metropolis-Hasting}
}

\section{Argument Details}{
  \emph{\code{Data  = list(X, share, J, Z)}} [\code{Z} optional]
  \tabular{ll}{
    \code{J:      } \tab number of alternatives, excluding an outside option \cr
    \code{X:      } \tab \eqn{J*T x K} matrix (no outside option, which is normalized to 0). \cr
    \code{        } \tab If IV is used, the last column of \code{X} is the endogeneous variable. \cr
    \code{share:  } \tab \eqn{J*T} vector (no outside option). \cr
    \code{        } \tab Note that both the \code{share} vector and the \code{X} matrix are organized by the \eqn{jt} index. \cr
    \code{        } \tab \eqn{j} varies faster than \eqn{t}, i.e. \eqn{(j=1,t=1), (j=2,t=1), ..., (j=J,T=1), ..., (j=J,t=T)} \cr
    \code{Z:      } \tab \eqn{J*T x I} matrix of instrumental variables (optional)
    }
  \emph{\code{Prior = list(sigmasqR, theta_hat, A, deltabar, Ad, nu0, s0_sq, VOmega)} [optional]}
  \tabular{ll}{
    \code{sigmasqR:   } \tab \eqn{K*(K+1)/2} vector for \eqn{r} prior variance (def: diffuse prior for \eqn{\Sigma}) \cr
    \code{theta_hat:  } \tab \eqn{K} vector for \eqn{\theta_bar} prior mean (def: 0 vector) \cr
    \code{A:          } \tab \eqn{K x K} matrix for \eqn{\theta_bar} prior precision (def: \code{0.01*diag(K)}) \cr
    \code{deltabar:   } \tab \eqn{I} vector for \eqn{\delta} prior mean (def: 0 vector) \cr
    \code{Ad:         } \tab \eqn{I x I} matrix for \eqn{\delta} prior precision (def: \code{0.01*diag(I)}) \cr
    \code{nu0:        } \tab d.f. parameter for \eqn{\tau_sq} and \eqn{\Omega} prior (def: K+1) \cr
    \code{s0_sq:      } \tab scale parameter for \eqn{\tau_sq} prior (def: 1) \cr
    \code{VOmega:     } \tab \eqn{2 x 2} matrix parameter for \eqn{\Omega} prior (def: \code{matrix(c(1,0.5,0.5,1),2,2)})
    }
  \emph{\code{Mcmc  = list(R, keep, nprint, H, initial_theta_bar, initial_r, initial_tau_sq, initial_Omega, initial_delta, s, cand_cov, tol)} [only \code{R} and \code{H} required]}
  \tabular{ll}{
    \code{R:                  } \tab number of MCMC draws \cr
    \code{keep:               } \tab MCMC thinning parameter -- keep every \code{keep}th draw (def: 1) \cr
    \code{nprint:             } \tab print the estimated time remaining for every \code{nprint}'th draw (def: 100, set to 0 for no print) \cr
    \code{H:                  } \tab number of random draws used for Monte-Carlo integration \cr
    \code{initial_theta_bar:  } \tab initial value of \eqn{\theta_bar} (def: 0 vector) \cr
    \code{initial_r:          } \tab initial value of \eqn{r} (def: 0 vector) \cr
    \code{initial_tau_sq:     } \tab initial value of \eqn{\tau_sq} (def: 0.1) \cr
    \code{initial_Omega:      } \tab initial value of \eqn{\Omega} (def: \code{diag(2)}) \cr
    \code{initial_delta:      } \tab initial value of \eqn{\delta} (def: 0 vector) \cr
    \code{s:                  } \tab scale parameter of Metropolis-Hasting increment (def: automatically tuned) \cr
    \code{cand_cov:           } \tab var-cov matrix of Metropolis-Hasting increment (def: automatically tuned) \cr
    \code{tol:                } \tab convergence tolerance for the contraction mapping (def: 1e-6)
    }
}

\section{Model Details}{
  \subsection{Model and Priors (without IV):}{ 
        \eqn{u_ijt = X_jt \theta_i + \eta_jt + e_ijt}\cr
        \eqn{e_ijt} \eqn{\sim}{~} type I Extreme Value (logit)\cr
        \eqn{\theta_i} \eqn{\sim}{~}  \eqn{N(\theta_bar, \Sigma)}\cr
        \eqn{\eta_jt} \eqn{\sim}{~} \eqn{N(0, \tau_sq)}\cr 
        
  This structure implies a logit model for each consumer (\eqn{\theta}). 
  Aggregate shares (\code{share}) are produced by integrating this consumer level 
  logit model over the assumed normal distribution of \eqn{\theta}. 

        \eqn{r} \eqn{\sim}{~} \eqn{N(0, diag(sigmasqR))}\cr
        \eqn{\theta_bar} \eqn{\sim}{~} \eqn{N(\theta_hat, A^-1)}\cr
        \eqn{\tau_sq} \eqn{\sim}{~} \eqn{nu0*s0_sq / \chi^2 (nu0)}\cr
        
        Note: we observe the aggregate level market share, not individual level choices.\cr
        
        Note: \eqn{r} is the vector of nonzero elements of cholesky root of \eqn{\Sigma}. 
        Instead of \eqn{\Sigma} we draw \eqn{r}, which is one-to-one correspondence with the positive-definite \eqn{\Sigma}.
  }
  \subsection{Model and Priors (with IV):}{ 
        \eqn{u_ijt = X_jt \theta_i + \eta_jt + e_ijt}\cr
        \eqn{e_ijt} \eqn{\sim}{~} type I Extreme Value (logit)\cr
        \eqn{\theta_i} \eqn{\sim}{~}  \eqn{N(\theta_bar, \Sigma)}\cr
        
        \eqn{X_jt = [X_exo_jt, X_endo_jt]}\cr
        \eqn{X_endo_jt = Z_jt \delta_jt + \zeta_jt}\cr
        \eqn{vec(\zeta_jt, \eta_jt)} \eqn{\sim}{~} \eqn{N(0, \Omega)}\cr
        
        \eqn{r} \eqn{\sim}{~} \eqn{N(0, diag(sigmasqR))}\cr
        \eqn{\theta_bar} \eqn{\sim}{~} \eqn{N(\theta_hat, A^-1)}\cr
        \eqn{\delta} \eqn{\sim}{~} \eqn{N(deltabar, Ad^-1)}\cr
        \eqn{\Omega} \eqn{\sim}{~} \eqn{IW(nu0, VOmega)}\cr
  }
}

\section{MCMC and Tuning Details:}{
   \subsection{MCMC Algorithm:}{
  Step 1 (\eqn{\Sigma}):\cr
        Given \eqn{\theta_bar} and \eqn{\tau_sq}, draw \eqn{r} via Metropolis-Hasting.\cr
        Covert the drawn \eqn{r} to \eqn{\Sigma}.\cr
        
        Note: if user does not specify the Metropolis-Hasting increment parameters 
        (\code{s} and \code{cand_cov}), \code{rbayesBLP} automatically tunes the parameters.
        
  Step 2 without IV (\eqn{\theta_bar}, \eqn{\tau_sq}):\cr
        Given \eqn{\Sigma}, draw \eqn{\theta_bar} and \eqn{\tau_sq} via Gibbs sampler.\cr
        
  Step 2 with IV (\eqn{\theta_bar}, \eqn{\delta}, \eqn{\Omega}):\cr
        Given \eqn{\Sigma}, draw \eqn{\theta_bar}, \eqn{\delta}, and \eqn{\Omega} via IV Gibbs sampler.\cr
  }
        
   \subsection{Tuning Metropolis-Hastings algorithm:}{
   r_cand = r_old + s*N(0,cand_cov)\cr
   Fix the candidate covariance matrix as cand_cov0 = diag(rep(0.1, K), rep(1, K*(K-1)/2)).\cr
   Start from s0 = 2.38/sqrt(dim(r))\cr
   
   Repeat\{\cr
   Run 500 MCMC chain.\cr   
   If acceptance rate < 30\% => update s1 = s0/5.\cr
   If acceptance rate > 50\% => update s1 = s0*3.\cr
   (Store r draws if acceptance rate is 20~80\%.)\cr
   s0 = s1\cr
   \} until acceptance rate is 30~50\%
   
   Scale matrix C = s1*sqrt(cand_cov0)\cr
   Correlation matrix R = Corr(r draws)\cr
   Use C*R*C as s^2*cand_cov.
   }
}

\references{For further discussion, see \emph{Bayesian Analysis of Random Coefficient Logit Models Using Aggregate Data} by Jiang, Manchanda, and Rossi, \emph{Journal of Econometrics}, 2009. \cr 
}

\author{Keunwoo Kim, Anderson School, UCLA, \email{keunwoo.kim@gmail.com}.}

\examples{
if(nchar(Sys.getenv("LONG_TEST")) != 0) {

## Simulate aggregate level data
simulData <- function(para, others, Hbatch) {
  # Hbatch does the integration for computing market shares
  #      in batches of size Hbatch

  ## parameters
  theta_bar <- para$theta_bar
  Sigma <- para$Sigma
  tau_sq <- para$tau_sq

  T <- others$T	
  J <- others$J	
  p <- others$p	
  H <- others$H	
  K <- J + p	

  ## build X	
  X <- matrix(runif(T*J*p), T*J, p)
  inter <- NULL
  for (t in 1:T) { inter <- rbind(inter, diag(J)) }
  X <- cbind(inter, X)

  ## draw eta ~ N(0, tau_sq)	
  eta <- rnorm(T*J)*sqrt(tau_sq)
  X <- cbind(X, eta)

  share <- rep(0, J*T)
  for (HH in 1:(H/Hbatch)){
    ## draw theta ~ N(theta_bar, Sigma)
    cho <- chol(Sigma)
    theta <- matrix(rnorm(K*Hbatch), nrow=K, ncol=Hbatch)
    theta <- t(cho)\%*\%theta + theta_bar

    ## utility
    V <- X\%*\%rbind(theta, 1)
    expV <- exp(V)
    expSum <- matrix(colSums(matrix(expV, J, T*Hbatch)), T, Hbatch)
    expSum <- expSum \%x\% matrix(1, J, 1)
    choiceProb <- expV / (1 + expSum)
    share <- share +  rowSums(choiceProb) / H
  }

  ## the last K+1'th column is eta, which is unobservable.
  X <- X[,c(1:K)]	
  return (list(X=X, share=share))
}

## true parameter
theta_bar_true <- c(-2, -3, -4, -5)
Sigma_true <- rbind(c(3,2,1.5,1), c(2,4,-1,1.5), c(1.5,-1,4,-0.5), c(1,1.5,-0.5,3))
cho <- chol(Sigma_true)
r_true <- c(log(diag(cho)), cho[1,2:4], cho[2,3:4], cho[3,4]) 
tau_sq_true <- 1

## simulate data
set.seed(66)
T <- 300
J <- 3
p <- 1
K <- 4
H <- 1000000
Hbatch <- 5000

dat <- simulData(para=list(theta_bar=theta_bar_true, Sigma=Sigma_true, tau_sq=tau_sq_true),
                 others=list(T=T, J=J, p=p, H=H), Hbatch)
X <- dat$X
share <- dat$share

## Mcmc run
R <- 2000
H <- 50
Data1 <- list(X=X, share=share, J=J)
Mcmc1 <- list(R=R, H=H, nprint=0)
set.seed(66)
out <- rbayesBLP(Data=Data1, Mcmc=Mcmc1)

## acceptance rate
out$acceptrate

## summary of draws
summary(out$thetabardraw)
summary(out$Sigmadraw)
summary(out$tausqdraw)

### plotting draws
plot(out$thetabardraw)
plot(out$Sigmadraw)
plot(out$tausqdraw)
}
}

\keyword{models}