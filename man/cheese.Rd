\name{cheese}
\alias{cheese}
\docType{data}
\title{ Sliced Cheese Data}
\description{
  Panel data with sales volume for a package of Borden Sliced Cheese
  as well as a measure of display activity and price.  Weekly data aggregated
  to the "key" account or retailer/market level.
}
\usage{data(cheese)}
\format{
  A data frame with 5555 observations on the following 4 variables.
  \describe{
    \item{\code{RETAILER}}{a list of 88 retailers}
    \item{\code{VOLUME}}{unit sales}
    \item{\code{DISP}}{a measure of display activity -- per cent ACV on display}
    \item{\code{PRICE}}{in \$}
  }
}
\source{
  Boatwright et al (1999), "Account-Level Modeling for Trade Promotion," 
  \emph{JASA} 94, 1063-1073.
}
\references{
 Chapter 3, \emph{Bayesian Statistics and Marketing} by Rossi et al. \cr
 \url{http://www.perossi.org/home/bsm-1l}
}
\examples{
data(cheese)
cat(" Quantiles of the Variables ",fill=TRUE)
mat=apply(as.matrix(cheese[,2:4]),2,quantile)
print(mat)

##
## example of processing for use with rhierLinearModel
##
if(0)
{

retailer=levels(cheese$RETAILER)
nreg=length(retailer)
nvar=3
regdata=NULL
for (reg in 1:nreg) {
	y=log(cheese$VOLUME[cheese$RETAILER==retailer[reg]])
	iota=c(rep(1,length(y)))
	X=cbind(iota,cheese$DISP[cheese$RETAILER==retailer[reg]],
		log(cheese$PRICE[cheese$RETAILER==retailer[reg]]))
	regdata[[reg]]=list(y=y,X=X)
}
Z=matrix(c(rep(1,nreg)),ncol=1)
nz=ncol(Z)
##
## run each individual regression and store results
##
lscoef=matrix(double(nreg*nvar),ncol=nvar)
for (reg in 1:nreg) {
	coef=lsfit(regdata[[reg]]$X,regdata[[reg]]$y,intercept=FALSE)$coef
	if (var(regdata[[reg]]$X[,2])==0)  { lscoef[reg,1]=coef[1]; lscoef[reg,3]=coef[2]}
	else {lscoef[reg,]=coef }
}

R=2000
Data=list(regdata=regdata,Z=Z)
Mcmc=list(R=R,keep=1)

set.seed(66)
out=rhierLinearModel(Data=Data,Mcmc=Mcmc)

cat("Summary of Delta Draws",fill=TRUE)
summary(out$Deltadraw)
cat("Summary of Vbeta Draws",fill=TRUE)
summary(out$Vbetadraw)

if(0){
#
# plot hier coefs
plot(out$betadraw)
}

}
}
\keyword{datasets}
