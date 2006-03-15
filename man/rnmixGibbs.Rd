\name{rnmixGibbs}
\alias{rnmixGibbs}
\concept{bayes}
\concept{MCMC}
\concept{normal mixtures}
\concept{Gibbs Sampling}

\title{ Gibbs Sampler for Normal Mixtures}
\description{
  \code{rnmixGibbs} implements a Gibbs Sampler for normal mixtures.
}
\usage{
rnmixGibbs(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y) }
  \item{Prior}{ list(Mubar,A,nu,V,a,ncomp) (only ncomp required)}
  \item{Mcmc}{ list(R,keep) (R required) }
}
\details{
  Model: \cr
        \eqn{y_i} \eqn{\sim}{~} \eqn{N(mu_{ind_i},Sigma_{ind_i})}. \cr
        ind \eqn{\sim}{~} iid multinomial(p).  p is a ncomp x 1 vector of probs. 

  Priors:\cr
        \eqn{mu_j} \eqn{\sim}{~} \eqn{N(mubar,Sigma_j (x) A^{-1})}. \eqn{mubar=vec(Mubar)}. \cr
        \eqn{Sigma_j} \eqn{\sim}{~} IW(nu,V).\cr
        note: this is the natural conjugate prior -- a special case of multivariate 
             regression.\cr
        \eqn{p} \eqn{\sim}{~} Dirchlet(a).

  Output of the components is in the form of a list of lists. \cr
  compsdraw[[i]] is ith draw -- list of ncomp lists. \cr
  compsdraw[[i]][[j]] is list of parms for jth normal component. \cr
  jcomp=compsdraw[[i]][j]]. Then jth comp \eqn{\sim}{~} \eqn{N(jcomp[[1]],Sigma)}, 
  \eqn{Sigma} = t(R)\%*\%R, \eqn{R^{-1}} = jcomp[[2]].

  List arguments contain:
  \itemize{
    \item{y}{ n x k array of data (rows are obs) }
    \item{Mubar}{ 1 x k array with prior mean of normal comp means (def: 0)}
    \item{A}{ 1 x 1 precision parameter for prior on mean of normal comp (def: .01)}
    \item{nu}{ d.f. parameter for prior on Sigma (normal comp cov matrix) (def: k+3)}
    \item{V}{ k x k location matrix of IW prior on Sigma (def: nuI)}
    \item{a}{ ncomp x 1 vector of Dirichlet prior parms (def: rep(5,ncomp))}
    \item{ncomp}{ number of normal components to be included }
    \item{R}{ number of MCMC draws }
    \item{keep}{ MCMC thinning parm: keep every keepth draw (def: 1)}
  }
}
\value{
  a list containing:
  \item{probdraw}{R/keep x ncomp array of mixture prob draws}
  \item{zdraw}{R/keep x nobs array of indicators of mixture comp identity for each obs}
  \item{compdraw}{R/keep lists of lists of comp parm draws}
}

\note{
  In this model, the component normal parameters are not-identified due to label-switching.
  However, the fitted mixture of normals density is identified as it is invariant to label-switching. 
  See Allenby et al, chapter 5 for details. Use \code{eMixMargDen} or \code{momMix} to compute
  posterior expectation or distribution of various identified parameters.
}

\references{ For further discussion, see \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch, Chapter 3. \cr
  \url{http://gsbwww.uchicago.edu/fac/peter.rossi/research/bsm.html}
}

\author{ Peter Rossi, Graduate School of Business, University of Chicago,
  \email{Peter.Rossi@ChicagoGsb.edu}.
}

\seealso{ \code{\link{rmixture}}, \code{\link{rmixGibbs}} ,\code{\link{eMixMargDen}}, \code{\link{momMix}},
 \code{\link{mixDen}}, \code{\link{mixDenBi}}}

\examples{
##
if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}

set.seed(66)
dim=5;  k=3   # dimension of simulated data and number of "true" components
sigma = matrix(rep(0.5,dim^2),nrow=dim);diag(sigma)=1
sigfac = c(1,1,1);mufac=c(1,2,3); compsmv=list()
for(i in 1:k) compsmv[[i]] = list(mu=mufac[i]*1:dim,sigma=sigfac[i]*sigma)
comps = list() # change to "rooti" scale
for(i in 1:k) comps[[i]] = list(mu=compsmv[[i]][[1]],rooti=solve(chol(compsmv[[i]][[2]])))
pvec=(1:k)/sum(1:k)

nobs=500
dm = rmixture(nobs,pvec,comps)

Data=list(y=dm$x)
ncomp=9
Prior=list(ncomp=ncomp)
Mcmc=list(R=R,keep=1)
out=rnmixGibbs(Data=Data,Prior=Prior,Mcmc=Mcmc)

tmom=momMix(matrix(pvec,nrow=1),list(comps))
if(R < 1000) {begin=1} else {begin=500}
pmom=momMix(out$probdraw[begin:R,],out$compdraw[begin:R])
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

if(0){
##
## plotting examples
##
## check true and estimated marginal densities
end=R
grid=NULL
for (i in 1:dim){
  rgi=range(dm$x[,i])
  gr=seq(from=rgi[1],to=rgi[2],length.out=50)
  grid=cbind(grid,gr)
}

tmden=mixDen(grid,pvec,comps)
pmden=eMixMargDen(grid,out$probdraw[begin:end,],out$compdraw[begin:end])

## plot the marginal on third variable
plot(range(grid[,3]),c(0,1.1*max(tmden[,3],pmden[,3])),type="n",xlab="",ylab="density")
lines(grid[,3],tmden[,3],col="blue",lwd=2)
lines(grid[,3],pmden[,3],col="red",lwd=2)

## compute implied bivariate marginal distributions
i=1
j=2
rxi=range(dm$x[,1])
rxj=range(dm$x[,2])
xi=seq(from=rxi[1],to=rxi[2],length.out=50)
xj=seq(from=rxj[1],to=rxj[2],length.out=50)
den=matrix(0,ncol=length(xi),nrow=length(xj))
for(ind in as.integer(seq(from=begin,to=end,length.out=100))){
  den=den+mixDenBi(i,j,xi,xj,out$probdraw[ind,],out$compdraw[[ind]])
}
tden=matrix(0,ncol=length(xi),nrow=length(xj))
tden=mixDenBi(i,j,xi,xj,pvec,comps)

par(mfrow=c(2,1))
image(xi,xj,tden,col=terrain.colors(100),xlab="",ylab="")
contour(xi,xj,tden,add=TRUE,drawlabels=FALSE)
title("True Bivariate Marginal")
image(xi,xj,den,col=terrain.colors(100),xlab="",ylab="")
contour(xi,xj,den,add=TRUE,drawlabels=FALSE)
title("Posterior Mean of Bivariate Marginal")
}

}
\keyword{ multivariate }
