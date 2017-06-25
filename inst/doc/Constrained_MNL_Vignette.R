## ---- echo = FALSE-------------------------------------------------------
library(bayesm)
knitr::opts_chunk$set(fig.align = "center",
                      fig.height = 3.5,
                      warning = FALSE,
                      error = FALSE,
                      message = FALSE)

## ------------------------------------------------------------------------
# define function
drawprior <- function (mubar_betak, nvar, ncomp, a, nu, Amu, V, ndraw) {
  betakstar <- double(ndraw)
  betak     <- double(ndraw)
  otherbeta <- double(ndraw)
  mubar     <- c(rep(0, nvar-1), mubar_betak)
  
  for(i in 1:ndraw) {
    comps=list()
    for(k in 1:ncomp) {
      Sigma <- rwishart(nu,chol2inv(chol(V)))$IW
      comps[[k]] <- list(mubar + t(chol(Sigma/Amu)) %*% rnorm(nvar), 
                         backsolve(chol(Sigma), diag(1,nvar)) )
    }
    pvec         <- rdirichlet(a)
    beta         <- rmixture(1,pvec,comps)$x
    betakstar[i] <- beta[nvar]
    betak[i]     <- -exp(beta[nvar])
    otherbeta[i] <- beta[1]
  }
  
  return(list(betakstar=betakstar, betak=betak, otherbeta=otherbeta))
}
set.seed(1234)

## ------------------------------------------------------------------------
# specify rhierMnlRwMixture defaults
mubar_betak <- 0
nvar  <- 10
ncomp <- 3
a     <- rep(5, ncomp)
nu    <- nvar + 3
Amu   <- 0.01
V     <- nu*diag(c(rep(1,nvar-1),1))
ndraw <- 10000
defaultprior <- drawprior(mubar_betak, nvar, ncomp, a, nu, Amu, V, ndraw)

## ---- fig.align='center', fig.height=3.5, results='hold'-----------------
# plot priors under defaults
par(mfrow=c(1,3))
trimhist <- -20
hist(defaultprior$betakstar, breaks=40, col="magenta", 
     main="Beta_k_star", xlab="", ylab="", yaxt="n")
hist(defaultprior$betak[defaultprior$betak>trimhist],
     breaks=40, col="magenta", main="Beta_k",
     xlab="", ylab="", yaxt="n", xlim=c(trimhist,0))
hist(defaultprior$otherbeta, breaks=40, col="magenta",
     main="Other Beta", xlab="", ylab="", yaxt="n")

## ------------------------------------------------------------------------
# adjust priors for constraints
mubar_betak <- 2
nvar  <- 10
ncomp <- 3 
a     <- rep(5, ncomp)
nu    <- nvar + 15
Amu   <- 0.1 
V     <- nu*diag(c(rep(4,nvar-1),0.1))
ndraw <- 10000
tightprior <- drawprior(mubar_betak, nvar, ncomp, a, nu, Amu, V, ndraw)

## ---- fig.align='center', fig.height=3.5, results='hold'-----------------
# plot priors under adjusted values
par(mfrow=c(1,3))
trimhist <- -20
hist(tightprior$betakstar, breaks=40, col="magenta", 
     main="Beta_k_star", xlab="", ylab="", yaxt="n")
hist(tightprior$betak[tightprior$betak>trimhist],
     breaks=40, col="magenta", main="Beta_k",
     xlab="", ylab="", yaxt="n", xlim=c(trimhist,0))
hist(tightprior$otherbeta, breaks=40, col="magenta", 
     main="Other Beta", xlab="", ylab="", yaxt="n")

## ------------------------------------------------------------------------
library(bayesm)
data(camera)
length(camera)
str(camera[[1]])

## ------------------------------------------------------------------------
str(camera[[1]]$y)
str(as.data.frame(camera[[1]]$X))

## ---- echo = FALSE-------------------------------------------------------
nvar  <- 10
mubar <- c(rep(0,nvar-1),2)
Amu   <- 0.1 
ncomp <- 5 
a     <- rep(5, ncomp)
nu    <- 25
V     <- nu*diag(c(rep(4,nvar-1),0.1))

## ---- eval = FALSE-------------------------------------------------------
#  SignRes <- c(rep(0,nvar-1),-1)
#  
#  data  <- list(lgtdata=camera, p=5)
#  prior <- list(mubar=mubar, Amu=Amu, ncomp=ncomp, a=a, nu=nu, V=V, SignRes=SignRes)
#  mcmc  <- list(R=1e4, nprint=0)
#  
#  out <- rhierMnlRwMixture(Data=data, Prior=prior, Mcmc=mcmc)

## ---- echo = FALSE-------------------------------------------------------
temp <- capture.output(
          {SignRes <- c(rep(0,nvar-1),-1);
           data  <- list(lgtdata = camera, p = 5);
           prior <- list(mubar=mubar, Amu=Amu, ncomp=ncomp, a=a, nu=nu, V=V, SignRes=SignRes);
           mcmc  <- list(R = 1e4, nprint = 0);
           out <- rhierMnlRwMixture(Data = data, Prior = prior, Mcmc = mcmc)}, 
        file = NULL)

## ------------------------------------------------------------------------
par(mfrow=c(1,3))
ind_hist <- function(mod, i) {
  hist(mod$betadraw[i , 10, ], breaks = seq(-14,0,0.5), 
     col = "dodgerblue3", border = "grey", yaxt = "n",
     xlim = c(-14,0), xlab = "", ylab = "", main = paste("Ind.",i))
}
ind_hist(out,1)
ind_hist(out,2)
ind_hist(out,3)

## ------------------------------------------------------------------------
par(mfrow=c(1,1))
hist(apply(out$betadraw[ , 10, ], 1, mean), 
     xlim = c(-20,0), breaks = 20, 
     col = "firebrick2", border = "gray", xlab = "", ylab = "",
     main = "Posterior Means for Ind. Price Params, 
             With Sign Constraint")


## ------------------------------------------------------------------------
data0  <- list(lgtdata = camera, p = 5)
prior0 <- list(ncomp = 5)
mcmc0  <- list(R = 1e4, nprint = 0)

## ---- eval = FALSE-------------------------------------------------------
#  out0 <- rhierMnlRwMixture(Data = data0, Prior = prior0, Mcmc = mcmc0)

## ---- echo = FALSE-------------------------------------------------------
temp <- capture.output(
          {out0 <- rhierMnlRwMixture(Data = data0, Prior = prior0, Mcmc = mcmc0)}, 
        file = NULL)

## ------------------------------------------------------------------------
par(mfrow=c(1,3))
ind_hist <- function(mod, i) {
  hist(mod$betadraw[i , 10, ], breaks = seq(-12,2,0.5), 
     col = "dodgerblue4", border = "grey", yaxt = "n",
     xlim = c(-12,2), xlab = "", ylab = "", main = paste("Ind.",i))
}
ind_hist(out0,1)
ind_hist(out0,2)
ind_hist(out0,3)

## ------------------------------------------------------------------------
par(mfrow=c(1,1))
hist(apply(out0$betadraw[ , 10, ], 1, mean), 
     xlim = c(-15,5), breaks = 20, 
     col = "firebrick3", border = "gray", 
     xlab = "", ylab = "",
     main = "Posterior Means for Ind. Price Params, 
             No Sign Constraint")

