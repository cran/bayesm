\name{tuna}
\alias{tuna}
\docType{data}
\title{Data on Canned Tuna Sales}
\description{
  Volume of canned tuna sales as well as a measure of display activity, log price and log wholesale price.  
  Weekly data aggregated to the chain level.  This data is extracted from the Dominick's Finer Foods database
  maintained by the Kilts center for marketing at the University of Chicago's Booth School of Business.
  Brands are seven of the top 10 UPCs in the canned tuna product category.
}
\usage{data(tuna)}
\format{
  A data frame with 338 observations on the following 30 variables.
  \describe{
  \item{\code{WEEK}}{a numeric vector}
    \item{\code{MOVE1}}{unit sales of Star Kist 6 oz.}
    \item{\code{MOVE2}}{unit sales of Chicken of the Sea 6 oz.}
    \item{\code{MOVE3}}{unit sales of Bumble Bee Solid 6.12 oz.}
    \item{\code{MOVE4}}{unit sales of Bumble Bee Chunk 6.12 oz.}
    \item{\code{MOVE5}}{unit sales of Geisha 6 oz.}
    \item{\code{MOVE6}}{unit sales of Bumble Bee Large Cans.}
    \item{\code{MOVE7}}{unit sales of HH Chunk Lite 6.5 oz.}
    \item{\code{NSALE1}}{a measure of display activity of Star Kist 6 oz.}
    \item{\code{NSALE2}}{a measure of display activity of Chicken of the Sea 6 oz.}
    \item{\code{NSALE3}}{a measure of display activity of Bumble Bee Solid 6.12 oz.}
    \item{\code{NSALE4}}{a measure of display activity of Bumble Bee Chunk 6.12 oz.}
    \item{\code{NSALE5}}{a measure of display activity of Geisha 6 oz.}
    \item{\code{NSALE6}}{a measure of display activity of Bumble Bee Large Cans.}
    \item{\code{NSALE7}}{a measure of display activity of HH Chunk Lite 6.5 oz.}
    \item{\code{LPRICE1}}{log of price of Star Kist 6 oz.}
    \item{\code{LPRICE2}}{log of price of Chicken of the Sea 6 oz.}
    \item{\code{LPRICE3}}{log of price of Bumble Bee Solid 6.12 oz.}
    \item{\code{LPRICE4}}{log of price of Bumble Bee Chunk 6.12 oz.}
    \item{\code{LPRICE5}}{log of price of Geisha 6 oz.}
    \item{\code{LPRICE6}}{log of price of Bumble Bee Large Cans.}
    \item{\code{LPRICE7}}{log of price of HH Chunk Lite 6.5 oz.}
    \item{\code{LWHPRIC1}}{log of wholesale price of Star Kist 6 oz.}
    \item{\code{LWHPRIC2}}{log of wholesale price of Chicken of the Sea 6 oz.}
    \item{\code{LWHPRIC3}}{log of wholesale price of Bumble Bee Solid 6.12 oz.}
    \item{\code{LWHPRIC4}}{log of wholesale price of Bumble Bee Chunk 6.12 oz.}
    \item{\code{LWHPRIC5}}{log of wholesale price of Geisha 6 oz.}
    \item{\code{LWHPRIC6}}{log of wholesale price of Bumble Bee Large Cans.}
    \item{\code{LWHPRIC7}}{log of wholesale price of HH Chunk Lite 6.5 oz.}
    \item{\code{FULLCUST}}{total customers visits}
  }
}
\source{
  Chevalier, A. Judith, Anil K. Kashyap and Peter E. Rossi
  (2003), "Why Don't Prices Rise During Periods of Peak Demand? Evidence from Scanner Data," 
  \emph{The American Economic Review} , 93(1), 15-37.
}
\references{
 Chapter 7, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby and McCulloch. \cr
 \url{http://www.perossi.org/home/bsm-1}
}
\examples{
data(tuna)
cat(" Quantiles of sales",fill=TRUE)
mat=apply(as.matrix(tuna[,2:5]),2,quantile)
print(mat)

##
## example of processing for use with rivGibbs
##
if(0)
{
  data(tuna)                          
  t = dim(tuna)[1]    
  customers = tuna[,30]                 
  sales = tuna[,2:8]                                                        
  lnprice = tuna[,16:22]      
  lnwhPrice= tuna[,23:29]      
  share=sales/mean(customers)
  shareout=as.vector(1-rowSums(share))
  lnprob=log(share/shareout)  

# create w matrix

  I1=as.matrix(rep(1, t))
  I0=as.matrix(rep(0, t))
  intercept=rep(I1, 4)
  brand1=rbind(I1, I0, I0, I0)
  brand2=rbind(I0, I1, I0, I0)
  brand3=rbind(I0, I0, I1, I0)
  w=cbind(intercept, brand1, brand2, brand3)  
  
## choose brand 1 to 4        
          
  y=as.vector(as.matrix(lnprob[,1:4])) 
  X=as.vector(as.matrix(lnprice[,1:4]))     
  lnwhPrice=as.vector(as.matrix (lnwhPrice[1:4]))   
  z=cbind(w, lnwhPrice)
                        
  Data=list(z=z, w=w, x=X, y=y)
  Mcmc=list(R=R, keep=1)
  set.seed(66)
  out=rivGibbs(Data=Data,Mcmc=Mcmc)

  cat(" betadraws ",fill=TRUE)
  summary(out$betadraw)


if(0){
## plotting examples
plot(out$betadraw)
}
}

  
}
\keyword{datasets}
