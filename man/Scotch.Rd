\name{Scotch}
\alias{Scotch}
\docType{data}
\title{ Survey Data on Brands of Scotch Consumed}
\description{
  from Simmons Survey.  Brands used in last year for those
  respondents who report consuming scotch. 
}
\usage{data(Scotch)}
\format{
  A data frame with 2218 observations on the following 21 variables.
  All variables are coded 1 if consumed in last year, 0 if not.
  \describe{
    \item{\code{Chivas.Regal}}{a numeric vector}
    \item{\code{Dewar.s.White.Label}}{a numeric vector}
    \item{\code{Johnnie.Walker.Black.Label}}{a numeric vector}
    \item{\code{J...B}}{a numeric vector}
    \item{\code{Johnnie.Walker.Red.Label}}{a numeric vector}
    \item{\code{Other.Brands}}{a numeric vector}
    \item{\code{Glenlivet}}{a numeric vector}
    \item{\code{Cutty.Sark}}{a numeric vector}
    \item{\code{Glenfiddich}}{a numeric vector}
    \item{\code{Pinch..Haig.}}{a numeric vector}
    \item{\code{Clan.MacGregor}}{a numeric vector}
    \item{\code{Ballantine}}{a numeric vector}
    \item{\code{Macallan}}{a numeric vector}
    \item{\code{Passport}}{a numeric vector}
    \item{\code{Black...White}}{a numeric vector}
    \item{\code{Scoresby.Rare}}{a numeric vector}
    \item{\code{Grants}}{a numeric vector}
    \item{\code{Ushers}}{a numeric vector}
    \item{\code{White.Horse}}{a numeric vector}
    \item{\code{Knockando}}{a numeric vector}
    \item{\code{the.Singleton}}{a numeric vector}
  }
}
\source{
  Edwards, Y. and G. Allenby (2003), "Multivariate Analysis of Multiple Response Data,"
   \emph{JMR}  40, 321-334.
}
\references{
  Chapter 4, \emph{Bayesian Statistics and Marketing} by Rossi et al.\cr 
  \url{http://www.perossi.org/home/bsm-1}
}
\examples{
data(Scotch)
cat(" Frequencies of Brands", fill=TRUE)
mat=apply(as.matrix(Scotch),2,mean)
print(mat)
##
## use Scotch data to run Multivariate Probit Model
##
if(0){
##

y=as.matrix(Scotch)
p=ncol(y); n=nrow(y)
dimnames(y)=NULL
y=as.vector(t(y))
y=as.integer(y)
I_p=diag(p)
X=rep(I_p,n)
X=matrix(X,nrow=p)
X=t(X)

R=2000
Data=list(p=p,X=X,y=y)
Mcmc=list(R=R)
set.seed(66)
out=rmvpGibbs(Data=Data,Mcmc=Mcmc)

ind=(0:(p-1))*p + (1:p)
cat(" Betadraws ",fill=TRUE)
mat=apply(out$betadraw/sqrt(out$sigmadraw[,ind]),2,quantile,probs=c(.01,.05,.5,.95,.99))
attributes(mat)$class="bayesm.mat"
summary(mat)
rdraw=matrix(double((R)*p*p),ncol=p*p)
rdraw=t(apply(out$sigmadraw,1,nmat))
attributes(rdraw)$class="bayesm.var"
cat(" Draws of Correlation Matrix ",fill=TRUE)
summary(rdraw)

}

}
\keyword{datasets}
