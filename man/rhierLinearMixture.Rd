\name{rhierLinearMixture}
\alias{rhierLinearMixture}
\concept{bayes}
\concept{MCMC}
\concept{Gibbs Sampling}
\concept{mixture of normals}
\concept{normal mixture}
\concept{heterogeneity}
\concept{regresssion}
\concept{hierarchical models}
\concept{linear model}

\title{ Gibbs Sampler for Hierarchical Linear Model }
\description{
  \code{rhierLinearMixture} implements a Gibbs Sampler for hierarchical linear models with a mixture of normals prior.
}
\usage{
rhierLinearMixture(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(regdata,Z) (Z optional). }
  \item{Prior}{ list(deltabar,Ad,mubar,Amu,nu,V,nu.e,ssq,ncomp)  (all but ncomp are optional).}
  \item{Mcmc}{ list(R,keep) (R required).}
}
\details{
  Model: length(regdata) regression equations. \cr
        \eqn{y_i = X_ibeta_i + e_i}. \eqn{e_i} \eqn{\sim}{~} \eqn{N(0,tau_i)}.  nvar X vars in each equation. 

 Priors:\cr
        \eqn{tau_i} \eqn{\sim}{~} nu.e*\eqn{ssq_i/\chi^2_{nu.e}}.  \eqn{tau_i} is the variance of \eqn{e_i}.\cr

        \eqn{beta_i}= ZDelta[i,] + \eqn{u_i}. \cr
        Note: here ZDelta refers to Z\%*\%D, ZDelta[i,] is ith row of this product.\cr
        Delta is an nz x nvar array. 

        \eqn{u_i} \eqn{\sim}{~} \eqn{N(mu_{ind},Sigma_{ind})}. \eqn{ind} \eqn{\sim}{~} multinomial(pvec). \cr

        \eqn{pvec} \eqn{\sim}{~} dirichlet (a)\cr
        \eqn{delta= vec(Delta)} \eqn{\sim}{~} \eqn{N(deltabar,A_d^{-1})}\cr
        \eqn{mu_j} \eqn{\sim}{~} \eqn{N(mubar,Sigma_j (x) Amu^{-1})}\cr
        \eqn{Sigma_j} \eqn{\sim}{~} IW(nu,V) \cr


  List arguments contain:
  \itemize{
    \item{\code{regdata}}{ list of lists with X,y matrices for each of length(regdata) regressions}
    \item{\code{regdata[[i]]$X}}{ X matrix for equation i }
    \item{\code{regdata[[i]]$y}}{ y vector for equation i }
    \item{\code{deltabar}}{nz*nvar vector of prior means (def: 0)}
    \item{\code{Ad}}{ prior prec matrix for vec(Delta) (def: .01I)}
    \item{\code{mubar}}{ nvar x 1 prior mean vector for normal comp mean (def: 0)}
    \item{\code{Amu}}{ prior precision for normal comp mean (def: .01I)}
    \item{\code{nu}}{ d.f. parm for IW prior on norm comp Sigma (def: nvar+3)}
    \item{\code{V}}{ pds location parm for IW prior on norm comp Sigma (def: nuI)}
    \item{\code{nu.e}}{ d.f. parm for regression error variance prior (def: 3)}
    \item{\code{ssq}}{ scale parm for regression error var prior (def: var(\eqn{y_i}))}
    \item{\code{ncomp}}{ number of components used in normal mixture }
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
   }
}
\value{
  a list containing
  \item{taudraw}{R/keep x nreg array of error variance draws}
  \item{betadraw}{nreg x nvar x R/keep array of individual regression coef draws}
  \item{Deltadraw}{R/keep x nz x nvar array of Deltadraws}
  \item{probdraw}{ R/keep x ncomp matrix of draws of probs of mixture components (pvec)}
  \item{compdraw}{ list of list of lists (length R/keep)}
}
\note{
  More on \code{compdraw} component of return value list: \cr
  \itemize{
  \item{compdraw[[i]]}{ the ith draw of components for mixtures}
  \item{compdraw[[i]][[j]]}{ ith draw of the jth normal mixture comp}
  \item{compdraw[[i]][[j]][[1]]}{ ith draw of jth normal mixture comp mean vector}
  \item{compdraw[[i]][[j]][[2]]}{ ith draw of jth normal mixture cov parm (rooti) }
  }

  Note: Z should \strong{not} include an intercept and should be centered for ease of interpretation.\cr
  
  Be careful in assessing prior parameter, Amu.  .01 can be too small for some applications. See 
  Rossi et al, chapter 5 for full discussion.\cr

} 
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 3. \cr
  \url{http://faculty.chicagogsb.edu/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}
\seealso{ \code{\link{rhierLinearModel}} }
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(66)
nreg=300; nobs=500; nvar=3; nz=2

Z=matrix(runif(nreg*nz),ncol=nz) 
Z=t(t(Z)-apply(Z,2,mean))
Delta=matrix(c(1,-1,2,0,1,0),ncol=nz)
tau0=.1
iota=c(rep(1,nobs))

## create arguments for rmixture

tcomps=NULL
a=matrix(c(1,0,0,0.5773503,1.1547005,0,-0.4082483,0.4082483,1.2247449),ncol=3)
tcomps[[1]]=list(mu=c(0,-1,-2),rooti=a) 
tcomps[[2]]=list(mu=c(0,-1,-2)*2,rooti=a)
tcomps[[3]]=list(mu=c(0,-1,-2)*4,rooti=a)
tpvec=c(.4,.2,.4)                               

regdata=NULL						  # simulated data with Z
betas=matrix(double(nreg*nvar),ncol=nvar)
tind=double(nreg)

# simulate datasets with/without Z, using same components and tau

for (reg in 1:nreg) {
tempout=rmixture(1,tpvec,tcomps)
betas[reg,]=Delta\%*\%Z[reg,]+as.vector(tempout$x)
tind[reg]=tempout$z
X=cbind(iota,matrix(runif(nobs*(nvar-1)),ncol=(nvar-1)))
tau=tau0*runif(1,min=0.5,max=1)
y=X\%*\%betas[reg,]+sqrt(tau)*rnorm(nobs)
regdata[[reg]]=list(y=y,X=X,beta=betas[reg,],tau=tau)
}

## run rhierLinearMixture

Data=list(regdata=regdata,Z=Z)
Prior=list(ncomp=3)
Mcmc=list(R=R,keep=1)

out1=rhierLinearMixture(Data=Data,Prior=Prior,Mcmc=Mcmc)

if(R>1000) {begin=501} else {begin=1}

apply(out1$Deltadraw[begin:R,],2,mean)
cat(" Deltadraws ",fill=TRUE)
mat=apply(out1$Deltadraw,2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(as.vector(Delta),mat)
rownames(mat)[1]="delta"
print(mat)

if(0){
## plotting examples 
## plot histograms of draws of betas for nth regression

betahist=function(n)
{
par(mfrow=c(1,3))
hist(out1$betadraw[n,1,begin:R],breaks=30,main="beta1 with Z",xlab=" ",ylab=" ")
abline(v=hhbetamean1[n,1],col="red",lwd=2)
abline(v=regdata[[n]]$beta[1],col="blue",lwd=2)
hist(out1$betadraw[n,2,begin:R],breaks=30,main="beta2 with Z",xlab=" ",ylab=" ")
abline(v=hhbetamean1[n,2],col="red",lwd=2)
abline(v=regdata[[n]]$beta[2],col="blue",lwd=2)
hist(out1$betadraw[n,3,begin:R],breaks=30,main="beta3 with Z",xlab=" ",ylab=" ")
abline(v=hhbetamean1[n,3],col="red",lwd=2)
abline(v=regdata[[n]]$beta[3],col="blue",lwd=2)
}


betahist(10)	## plot betas for 10th regression, using regdata
betahist(20)
betahist(30)

## plot univariate marginal density of betas

grid=NULL
for (i in 1:nvar){
  rgi=range(betas[,i])
  gr=seq(from=rgi[1],to=rgi[2],length.out=50)
  grid=cbind(grid,gr)
}

tmden=mixDen(grid,tpvec,tcomps)
pmden=eMixMargDen(grid,as.matrix(out1$probdraw[begin:R,]),out1$compdraw[begin:R])

par(mfrow=c(1,3))

for (i in 1:nvar){
plot(range(grid[,i]),c(0,1.1*max(tmden[,i],pmden[,i])),type="n",xlab="",ylab="density")
lines(grid[,i],tmden[,i],col="blue",lwd=2)
lines(grid[,i],pmden[,i],col="red",lwd=2)
}


# plot bivariate marginal density of betas

end=R
rx1=range(betas[,1])
rx2=range(betas[,2])
rx3=range(betas[,3])
x1=seq(from=rx1[1],to=rx1[2],length.out=50)
x2=seq(from=rx2[1],to=rx2[2],length.out=50)
x3=seq(from=rx3[1],to=rx3[2],length.out=50)
den12=matrix(0,ncol=length(x1),nrow=length(x2))
den23=matrix(0,ncol=length(x2),nrow=length(x3))
den13=matrix(0,ncol=length(x1),nrow=length(x3))

for(ind in as.integer(seq(from=begin,to=end,length.out=100))){
den12=den12+mixDenBi(1,2,x1,x2,as.matrix(out1$probdraw[ind,]),out1$compdraw[[ind]])
den23=den23+mixDenBi(2,3,x2,x3,as.matrix(out1$probdraw[ind,]),out1$compdraw[[ind]])
den13=den13+mixDenBi(1,3,x1,x3,as.matrix(out1$probdraw[ind,]),out1$compdraw[[ind]])
}

tden12=matrix(0,ncol=length(x1),nrow=length(x2))
tden23=matrix(0,ncol=length(x2),nrow=length(x3))
tden13=matrix(0,ncol=length(x1),nrow=length(x3))
tden12=mixDenBi(1,2,x1,x2,tpvec,tcomps)
tden23=mixDenBi(2,3,x2,x3,tpvec,tcomps)
tden13=mixDenBi(1,3,x1,x3,tpvec,tcomps)

par(mfrow=c(2,3))
image(x1,x2,tden12,col=terrain.colors(100),xlab="",ylab="")
contour(x1,x2,tden12,add=TRUE,drawlabels=FALSE)
title("True vars 1&2")
image(x2,x3,tden23,col=terrain.colors(100),xlab="",ylab="")
contour(x2,x3,tden23,add=TRUE,drawlabels=FALSE)
title("True vars 2&3")
image(x1,x3,tden13,col=terrain.colors(100),xlab="",ylab="")
contour(x1,x3,tden13,add=TRUE,drawlabels=FALSE)
title("True vars 1&3")
image(x1,x2,den12,col=terrain.colors(100),xlab="",ylab="")
contour(x1,x2,den12,add=TRUE,drawlabels=FALSE)
title("Posterior vars 1&2")
image(x2,x3,den23,col=terrain.colors(100),xlab="",ylab="")
contour(x2,x3,den23,add=TRUE,drawlabels=FALSE)
title("Posterior vars 2&3")
image(x1,x3,den13,col=terrain.colors(100),xlab="",ylab="")
contour(x1,x3,den13,add=TRUE,drawlabels=FALSE)
title("Posterior vars 1&3")
}
}
\keyword{ regression }
