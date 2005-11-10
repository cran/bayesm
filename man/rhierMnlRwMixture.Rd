\name{rhierMnlRwMixture}
\alias{rhierMnlRwMixture}
\concept{bayes}
\concept{MCMC}
\concept{Multinomial Logit}
\concept{mixture of normals}
\concept{normal mixture}
\concept{heterogeneity}
\concept{hierarchical models}

\title{ MCMC Algorithm for Hierarchical Multinomial Logit with Mixture of Normals Heterogeneity}
\description{
  \code{rhierMnlRwMixture} is a MCMC algorithm for a hierarchical multinomial logit with a mixture of normals 
  heterogeneity distribution.  This is a hybrid Gibbs Sampler with a RW Metropolis step for the MNL 
  coefficients for each panel unit.
}
\usage{
rhierMnlRwMixture(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(p,lgtdata,Z) ( Z is optional) }
  \item{Prior}{ list(deltabar,Ad,mubar,Amu,nu,V,ncomp) (all but ncomp are optional)}
  \item{Mcmc}{ list(s,w,R,keep) (R required)}
}
\details{
  Model: \cr
  \eqn{y_i} \eqn{\sim}{~} \eqn{MNL(X_i,theta_i)}.  i=1,\ldots, length(lgtdata). \eqn{theta_i} is nvar x 1.

  \eqn{theta_i}= ZDelta[i,] + \eqn{u_i}. \cr
  Note: here ZDelta refers to Z\%*\%D, ZDelta[i,] is ith row of this product.\cr
  Delta is an nz x nvar array. 

  \eqn{u_i} \eqn{\sim}{~} \eqn{N(mu_{ind},Sigma_{ind})}. \eqn{ind} \eqn{\sim}{~} multinomial(pvec). \cr

  Priors: \cr
  \eqn{pvec} \eqn{\sim}{~} dirichlet (a)\cr
  \eqn{delta= vec(Delta)} \eqn{\sim}{~} \eqn{N(deltabar,A_d^{-1})}\cr
  \eqn{mu_j} \eqn{\sim}{~} \eqn{N(mubar,Sigma_j (x) Amu^{-1})}\cr
  \eqn{Sigma_j} \eqn{\sim}{~} IW(nu,V) \cr

  Lists contain:
  \itemize{
    \item{\code{p}}{ p is number of choice alternatives}
    \item{\code{lgtdata}}{list of lists with each cross-section unit MNL data}
    \item{\code{lgtdata[[i]]$y}}{ \eqn{n_i} vector of multinomial outcomes (1,\ldots,m)}
    \item{\code{lgtdata[[i]]$X}}{ \eqn{n_i} by nvar design matrix for ith unit}
    \item{\code{deltabar}}{nz*nvar vector of prior means (def: 0)}
    \item{\code{Ad}}{ prior prec matrix for vec(D) (def: .01I)}
    \item{\code{mubar}}{ nvar x 1 prior mean vector for normal comp mean (def: 0)}
    \item{\code{Amu}}{ prior precision for normal comp mean (def: .01I)}
    \item{\code{nu}}{ d.f. parm for IW prior on norm comp Sigma (def: nvar+3)}
    \item{\code{V}}{ pds location parm for IW prior on norm comp Sigma (def: nuI)}
    \item{\code{ncomp}}{ number of components used in normal mixture }
    \item{\code{s}}{ scaling parm for RW Metropolis (def: 2.93/sqrt(nvar))}
    \item{\code{w}}{ fractional likelihood weighting parm (def: .1)}
    \item{\code{R}}{ number of MCMC draws}
    \item{\code{keep}}{ MCMC thinning parm: keep every keepth draw (def: 1)}
  }
}
\value{
  a list containing:
  \item{Deltadraw}{R/keep  x nz*nvar matrix of draws of Delta, first row is initial value}
  \item{betadraw}{ nlgt x nvar x R/keep array of draws of betas}
  \item{probdraw}{ R/keep x ncomp matrix of draws of probs of mixture components (pvec)}
  \item{compdraw}{ list of list of lists (length R/keep)}
  \item{loglike}{ log-likelihood for each kept draw (length R/keep)}
}
\note{
  More on \code{compdraw} component of return value list: \cr
  \itemize{
  \item{compdraw[[i]]}{ the ith draw of components for mixtures}
  \item{compdraw[[i]][[j]]}{ ith draw of the jth normal mixture comp}
  \item{compdraw[[i]][[j]][[1]]}{ ith draw of jth normal mixture comp mean vector}
  \item{compdraw[[i]][[j]][[2]]}{ ith draw of jth normal mixture cov parm (rooti) }
  }

  Note: Z does \strong{not} include an intercept and is centered for ease of interpretation.\cr
  
  Be careful in assessing prior parameter, Amu.  .01 is too small for many applications. See 
  Rossi et al, chapter 5 for full discussion.\cr

  Note: as of version 2.0-2 of \code{bayesm}, the fractional weight parameter has been changed
  to a weight between 0 and 1.  w is the fractional weight on the normalized pooled likelihood.
  This differs from what is in Rossi et al chapter 5, i.e.

  \eqn{like_i^(1-w) x like_pooled^(n_i/N)*w}

 
  Large R values may be required (>20,000).

} 
\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 5. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}
  
\seealso{ \code{\link{rmnlIndepMetrop}} }
\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) 
{

set.seed(66)
p=3                                # num of choice alterns
ncoef=3  
nlgt=300                           # num of cross sectional units
nz=2
Z=matrix(runif(nz*nlgt),ncol=nz)
Z=t(t(Z)-apply(Z,2,mean))          # demean Z
ncomp=3                                # no of mixture components
Delta=matrix(c(1,0,1,0,1,2),ncol=2)
comps=NULL
comps[[1]]=list(mu=c(0,-1,-2),rooti=diag(rep(1,3)))
comps[[2]]=list(mu=c(0,-1,-2)*2,rooti=diag(rep(1,3)))
comps[[3]]=list(mu=c(0,-1,-2)*4,rooti=diag(rep(1,3)))
pvec=c(.4,.2,.4)

## simulate data
simlgtdata=NULL
ni=rep(50,300)
for (i in 1:nlgt) 
{  betai=Delta\%*\%Z[i,]+as.vector(rmixture(1,pvec,comps)$x)
   X=NULL
   for(j in 1:ni[i]) 
      { Xone=cbind(rbind(c(rep(0,p-1)),diag(p-1)),runif(p,min=-1.5,max=0))
        X=rbind(X,Xone) }
   outa=simmnlwX(ni[i],X,betai)
   simlgtdata[[i]]=list(y=outa$y,X=X,beta=betai)
}

## plot betas
if(0){
## set if(1) above to produce plots
bmat=matrix(0,nlgt,ncoef)
for(i in 1:nlgt) {bmat[i,]=simlgtdata[[i]]$beta}
par(mfrow=c(ncoef,1))
for(i in 1:ncoef) hist(bmat[,i],breaks=30,col="magenta")
}

##   set parms for priors and Z
nu=ncoef+3
V=nu*diag(rep(1,ncoef))
Ad=.01*(diag(rep(1,nz*ncoef)))
mubar=matrix(rep(0,ncoef),nrow=1)
deltabar=rep(0,ncoef*nz)
Amu=matrix(.01,ncol=1)
a=rep(5,ncoef)

R=10000
keep=5
w=.1
s=2.93/sqrt(ncoef)
Data1=list(p=p,lgtdata=simlgtdata,Z=Z)
Prior1=list(ncomp=ncomp,nu=nu,V=V,Amu=Amu,mubar=mubar,a=a,Ad=Ad,deltabar=deltabar)
Mcmc1=list(s=s,w=w,R=R,keep=keep)
out=rhierMnlRwMixture(Data=Data1,Prior=Prior1,Mcmc=Mcmc1)

if(R < 1000) {begin=1} else {begin=1000/keep}
end=R/keep
cat(" Deltadraws ",fill=TRUE)
mat=apply(out$Deltadraw[begin:end,],2,quantile,probs=c(.01,.05,.5,.95,.99))
mat=rbind(as.vector(Delta),mat); rownames(mat)[1]="delta"; print(mat)

tmom=momMix(matrix(pvec,nrow=1),list(comps))
pmom=momMix(out$probdraw[begin:end,],out$compdraw[begin:end])
mat=rbind(tmom$mu,pmom$mu)
rownames(mat)=c("true","post expect")
cat(" mu and posterior expectation of mu",fill=TRUE)
print(mat)
mat=rbind(tmom$sd,pmom$sd)
rownames(mat)=c("true","post expect")
cat(" std dev and posterior expectation of sd",fill=TRUE)
print(mat)
mat=rbind(as.vector(tmom$corr),as.vector(pmom$corr))
rownames(mat)=c("true","post expect")
cat(" corr and posterior expectation of corr",fill=TRUE)
print(t(mat)) 

if(0) {
## set if(1) to produce plots
par(mfrow=c(4,1))
plot(out$betadraw[1,1,])
abline(h=simlgtdata[[1]]$beta[1])
plot(out$betadraw[2,1,])
abline(h=simlgtdata[[2]]$beta[1])
plot(out$betadraw[100,3,])
abline(h=simlgtdata[[100]]$beta[3])
plot(out$betadraw[101,3,])
abline(h=simlgtdata[[101]]$beta[3])
par(mfrow=c(4,1))
plot(out$Deltadraw[,1])
abline(h=Delta[1,1])
plot(out$Deltadraw[,2])
abline(h=Delta[2,1])
plot(out$Deltadraw[,3])
abline(h=Delta[3,1])
plot(out$Deltadraw[,6])
abline(h=Delta[3,2])
begin=1000/keep
end=R/keep
ngrid=50
grid=matrix(0,ngrid,ncoef)
rgm=matrix(c(-3,-7,-10,3,1,0),ncol=2)
for(i in 1:ncoef) {rg=rgm[i,]; grid[,i]=rg[1] + ((1:ngrid)/ngrid)*(rg[2]-rg[1])}
mden=eMixMargDen(grid,out$probdraw[begin:end,],out$compdraw[begin:end])
par(mfrow=c(2,ncoef))
for(i in 1:ncoef)
{plot(grid[,i],mden[,i],type="l")}
for(i in 1:ncoef)
tden=mixDen(grid,pvec,comps)
for(i in 1:ncoef)
{plot(grid[,i],tden[,i],type="l")}
}
}

}
\keyword{ models }
