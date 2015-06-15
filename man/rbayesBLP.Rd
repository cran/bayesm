\name{rbayesBLP}
\alias{rbayesBLP}
\concept{bayes}
\concept{random coefficient logit}
\concept{BLP}
\concept{Metropolis Hasting}

\title{ Bayesian Analysis of Random Coefficient Logit Models Using Aggregate Data }
\description{
  \code{rbayesBLP} implements a hybrid MCMC algorithm for aggregate level sales data in a market with differentiated products. Version 3.0-1 contains an error for use of instruments with this function. This will be fixed in version 3.0-2.
}
\usage{
rbayesBLP(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(X,share,J,Z) (X, share, and J: required). }
  \item{Prior}{ list(sigmasqR,theta_hat,A,deltabar,Ad,nu0,s0_sq,VOmega) (optional).}
  \item{Mcmc}{ list(R,H,initial_theta_bar,initial_r,initial_tau_sq,initial_Omega,initial_delta,s,cand_cov,tol,keep,nprint) (R and H: required).}
}
\details{
  Model: \cr
        \eqn{u_ijt = X_jt \theta_i + \eta_jt + e_ijt}\cr
        \eqn{e_ijt} \eqn{\sim}{~} type I Extreme Value (logit)\cr
        \eqn{\theta_i} \eqn{\sim}{~}  \eqn{N(\theta_bar, \Sigma)}\cr
        \eqn{\eta_jt} \eqn{\sim}{~} \eqn{N(0, \tau_sq)}\cr 
        This structure implies a logit model for each consumer (\eqn{\theta}). Aggregate shares \code{share} are produced by integrating this consumer level logit model over the assumed normal distribution of \eqn{\theta}. 

 Priors:\cr
        \eqn{r} \eqn{\sim}{~} \eqn{N(0,diag(sigmasqR))}.\cr
        \eqn{\theta_bar} \eqn{\sim}{~} \eqn{N(\theta_hat,A^-1)}.\cr
        \eqn{\tau_sq} \eqn{\sim}{~} \eqn{nu0*s0_sq / \chi^2 (nu0)}\cr
        
        Note: we observe the aggregate level market share, not individual level choice.\cr
        
        Note: \eqn{r} is the vector of nonzero elements of cholesky root of \eqn{\Sigma}. Instead of \eqn{\Sigma} we draw \eqn{r}, which is one-to-one correspondence with the positive-definite \eqn{\Sigma}.
        
 Model (with IV): \cr
        \eqn{u_ijt = X_jt \theta_i + \eta_jt + e_ijt}\cr
        \eqn{e_ijt} \eqn{\sim}{~} type I Extreme Value (logit)\cr
        \eqn{\theta_i} \eqn{\sim}{~}  \eqn{N(\theta_bar, \Sigma)}\cr
        
        \eqn{X_jt = [X_exo_jt, X_endo_jt]}\cr
        \eqn{X_endo_jt = Z_jt \delta_jt + \zeta_jt}\cr
        \eqn{vec(\zeta_jt, \eta_jt)} \eqn{\sim}{~} \eqn{N(0, \Omega)}\cr

 Priors (with IV):\cr
        \eqn{r} \eqn{\sim}{~} \eqn{N(0,diag(sigmasqR))}.\cr
        \eqn{\theta_bar} \eqn{\sim}{~} \eqn{N(\theta_hat,A^-1)}.\cr
        \eqn{\delta} \eqn{\sim}{~} \eqn{N(deltabar,Ad^-1)}.\cr
        \eqn{\Omega} \eqn{\sim}{~} \eqn{IW(nu0, VOmega)}\cr
          
 Step 1 (\eqn{\Sigma}):\cr
        Given \eqn{\theta_bar} and \eqn{\tau_sq}, draw \eqn{r} via Metropolis-Hasting.\cr
        Covert the drawn \eqn{r} to \eqn{\Sigma}.\cr
        
        Note: if user does not specify the Metropolis-Hasting increment parameters (\code{s} and \code{cand_cov}), \code{rbayesBLP} automatically tunes the parameters.
        
 Step 2 (\eqn{\theta_bar}, \eqn{\tau_sq}):\cr
        Given \eqn{\Sigma}, draw \eqn{\theta_bar} and \eqn{\tau_sq} via Gibbs sampler.\cr
        
 Step 2 (with IV: \eqn{\theta_bar}, \eqn{\delta}, \eqn{\Omega}):\cr
        Given \eqn{\Sigma}, draw \eqn{\theta_bar}, \eqn{\delta}, and \eqn{\Omega} via IV Gibbs sampler.\cr
        
  List arguments contain:\cr
  
  Data
  \itemize{
    \item{\code{J}}{ number of alternatives without outside option}
    \item{\code{X}}{ J*T by K matrix (no outside option, which is normalized to 0). If IV is used, the last column is endogeneous variable.}
    \item{\code{share}}{ J*T vector (no outside option)}
    \item{\code{Z}}{ J*T by I matrix of instrumental variables (optional)}
  }
  Note: both the \code{share} vector and the \code{X} matrix are organized by the jt index. j varies faster than t, i.e. (j=1,t=1),(j=2,t=1), ..., (j=J,T=1), ..., (j=J,t=T)\cr
  
  Prior
  \itemize{
    \item{\code{sigmasqR}}{ K*(K+1)/2 vector for \eqn{r} prior variance (def: diffuse prior for \eqn{\Sigma})}    
    \item{\code{theta_hat}}{ K vector for \eqn{\theta_bar} prior mean (def: 0 vector)}
    \item{\code{A}}{ K by K matrix for \eqn{\theta_bar} prior precision (def: 0.01*diag(K))}
    \item{\code{deltabar}}{ I vector for \eqn{\delta} prior mean (def: 0 vector)}
    \item{\code{Ad}}{ I by I matrix for \eqn{\delta} prior precision (def: 0.01*diag(I))}
    \item{\code{nu0}}{ d.f. parameter for \eqn{\tau_sq} and \eqn{\Omega} prior (def: K+1)}
    \item{\code{s0_sq}}{ scale parameter for \eqn{\tau_sq} prior (def: 1)}
    \item{\code{VOmega}}{ 2 by 2 matrix parameter for \eqn{\Omega} prior (def: matrix(c(1,0.5,0.5,1),2,2))}
  }
  Mcmc
  \itemize{    
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{H}}{ number of random draws used for Monte-Carlo integration}
    \item{\code{initial_theta_bar}}{ initial value of \eqn{\theta_bar} (def: 0 vector)}
    \item{\code{initial_r}}{ initial value of \eqn{r} (def: 0 vector)}
    \item{\code{initial_tau_sq}}{ initial value of \eqn{\tau_sq} (def: 0.1)}
    \item{\code{initial_Omega}}{ initial value of \eqn{\Omega} (def: diag(2))}
    \item{\code{initial_delta}}{ initial value of \eqn{\delta} (def: 0 vector)}
    \item{\code{s}}{ scale parameter of Metropolis-Hasting increment (def: automatically tuned)}
    \item{\code{cand_cov}}{ var-cov matrix of Metropolis-Hasting increment (def: automatically tuned)}
    \item{\code{tol}}{ convergence tolerance for the contraction mapping (def: 1e-6)}
    \item{\code{keep}}{ MCMC thinning parameter: keep every keepth draw (def: 1)}
    \item{\code{nprint}}{ print the estimated time remaining for every nprint'th draw (def: 100)}    
   }
   Tuning Metropolis-Hastings algorithm:\cr
   
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
   Use C*R*C as s^2*cand_cov.\cr
   
   
   
}
\value{
  a list containing
  \item{thetabardraw}{K by R/keep matrix of random coefficient mean draws}
  \item{Sigmadraw}{K*K by R/keep matrix of random coefficient variance draws}
  \item{rdraw}{K*K by R/keep matrix of \eqn{r} draws (same information as in \code{Sigmadraw})}
  \item{tausqdraw}{R/keep vector of aggregate demand shock variance draws}
  \item{Omegadraw}{2*2 by R/keep matrix of correlated endogenous shock variance draws}
  \item{deltadraw}{I by R/keep matrix of endogenous structural equation coefficient draws}
  \item{acceptrate}{scalor of acceptance rate of Metropolis-Hasting}
  \item{s}{scale parameter used for Metropolis-Hasting}
  \item{cand_cov}{var-cov matrix used for Metropolis-Hasting}
}
\references{ For further discussion, see \emph{Bayesian Analysis of Random Coefficient Logit Models Using Aggregate Data}
  by Jiang, Manchanda and Rossi, Journal of Econometrics, 2009. \cr
  \url{http://www.sciencedirect.com/science/article/pii/S0304407608002297}
}

\author{ Keunwoo Kim, Anderson School, UCLA,
  \email{keunwoo.kim@gmail.com}.
}

\examples{

if(nchar(Sys.getenv("LONG_TEST")) != 0) {
###
### Simulate aggregate level data
###
simulData <- function(para, others, Hbatch){
  #
  # Keunwoo Kim, UCLA Anderson
  #
  ### parameters
  theta_bar <- para$theta_bar
  Sigma <- para$Sigma
  tau_sq <- para$tau_sq
	
  T <- others$T	
  J <- others$J	
  p <- others$p	
  H <- others$H	
  K <- J + p	
  
  # Hbatch does the integration for computing market shares in batches of
  #        size Hbatch

  ### build X	
  X <- matrix(runif(T*J*p), T*J, p)
  inter <- NULL
  for (t in 1:T){
    inter <- rbind(inter, diag(J))
  }
  X <- cbind(inter, X)

  ### draw eta ~ N(0, tau_sq)	
  eta <- rnorm(T*J)*sqrt(tau_sq)
  X <- cbind(X, eta)
	
  share <- rep(0, J*T)
  for (HH in 1:(H/Hbatch)){
    ### draw theta ~ N(theta_bar, Sigma)
    cho <- chol(Sigma)
    theta <- matrix(rnorm(K*Hbatch), nrow=K, ncol=Hbatch)
    theta <- t(cho)\%*\%theta + theta_bar

    ### utility
    V <- X\%*\%rbind(theta, 1)
    expV <- exp(V)
    expSum <- matrix(colSums(matrix(expV, J, T*Hbatch)), T, Hbatch)
    expSum <- expSum \%x\% matrix(1, J, 1)
    choiceProb <- expV / (1 + expSum)
    share <- share +  rowSums(choiceProb) / H
  }	
	
  ### the last K+1'th column is eta, which is unobservable.
  X	<- X[,c(1:K)]	
  return (list(X=X, share=share))
}

### true parameter
theta_bar_true <- c(-2, -3, -4, -5)
Sigma_true <- rbind(c(3,2,1.5,1),c(2,4,-1,1.5),c(1.5,-1,4,-0.5),c(1,1.5,-0.5,3))
cho <- chol(Sigma_true)
r_true <- c(log(diag(cho)),cho[1,2:4],cho[2,3:4],cho[3,4]) 
tau_sq_true <- 1

### simulate data
set.seed(66)
T <- 300;J <- 3;p <- 1;K <- 4;H <- 1000000;Hbatch <- 5000
dat <- simulData(para=list(theta_bar=theta_bar_true, Sigma=Sigma_true, tau_sq=tau_sq_true),
        others=list(T=T, J=J, p=p, H=H), Hbatch)
X <- dat$X
share <- dat$share

### Mcmc run
R <- 2000;H <- 50
Data1 <- list(X=X, share=share, J=J)
Mcmc1 <- list(R=R, H=H, nprint=0)
set.seed(66)
out <- rbayesBLP(Data=Data1, Mcmc=Mcmc1)

### acceptance rate
out$acceptrate

### summary of draws
summary(out$thetabardraw)
summary(out$Sigmadraw)
summary(out$tausqdraw)

### plotting draws
plot(out$thetabardraw)
plot(out$Sigmadraw)
plot(out$tausqdraw)
}
}

