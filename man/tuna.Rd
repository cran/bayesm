\name{tuna}
\alias{tuna}
\docType{data}

\title{Canned Tuna Sales Data}

\description{
Volume of canned tuna sales as well as a measure of display activity, log price, and log wholesale price. Weekly data aggregated to the chain level. This data is extracted from the Dominick's Finer Foods database maintained by the Kilts Center for Marketing at the University of Chicago's Booth School of Business. Brands are seven of the top 10 UPCs in the canned tuna product category.
}

\usage{data(tuna)}

\format{
A data frame with 338 observations on 30 variables.
  \tabular{ll}{
    \ldots\code{$WEEK     } \tab a numeric vector \cr
    \ldots\code{$MOVE#    } \tab unit sales of brand # \cr
    \ldots\code{$NSALE#   } \tab a measure of display activity of brand # \cr
    \ldots\code{$LPRICE#  } \tab log of price of brand # \cr
    \ldots\code{$LWHPRIC# } \tab log of wholesale price of brand # \cr
    \ldots\code{$FULLCUST } \tab total customers visits 
  }
  
The brands are:
\tabular{ll}{
    1. \tab Star Kist 6 oz. \cr
    2. \tab Chicken of the Sea 6 oz. \cr
    3. \tab Bumble Bee Solid 6.12 oz. \cr
    4. \tab Bumble Bee Chunk 6.12 oz. \cr
    5. \tab Geisha 6 oz. \cr
    6. \tab Bumble Bee Large Cans. \cr
    7. \tab HH Chunk Lite 6.5 oz. 
  }
}

\source{Chevalier, Judith, Anil Kashyap, and Peter Rossi (2003), "Why Don't Prices Rise During Periods of Peak Demand? Evidence from Scanner Data," \emph{The American Economic Review} , 93(1), 15--37.}

\references{Chapter 7, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
data(tuna)
cat(" Quantiles of sales", fill=TRUE)
mat = apply(as.matrix(tuna[,2:5]), 2, quantile)
print(mat)


## example of processing for use with rivGibbs
if(0) {
  data(tuna)                          
  t = dim(tuna)[1]    
  customers = tuna[,30]                 
  sales = tuna[,2:8]                                                        
  lnprice = tuna[,16:22]      
  lnwhPrice = tuna[,23:29]      
  share = sales/mean(customers)
  shareout = as.vector(1-rowSums(share))
  lnprob = log(share/shareout)  

  ## create w matrix
  I1 = as.matrix(rep(1,t))
  I0 = as.matrix(rep(0,t))
  intercept = rep(I1,4)
  brand1 = rbind(I1, I0, I0, I0)
  brand2 = rbind(I0, I1, I0, I0)
  brand3 = rbind(I0, I0, I1, I0)
  w = cbind(intercept, brand1, brand2, brand3)  
  
  ## choose brand 1 to 4
  y = as.vector(as.matrix(lnprob[,1:4]))
  X = as.vector(as.matrix(lnprice[,1:4]))
  lnwhPrice = as.vector(as.matrix(lnwhPrice[1:4]))
  z = cbind(w, lnwhPrice)
                        
  Data = list(z=z, w=w, x=X, y=y)
  Mcmc = list(R=R, keep=1)
  
  set.seed(66)
  out = rivGibbs(Data=Data, Mcmc=Mcmc)

  cat(" betadraws ", fill=TRUE)
  summary(out$betadraw)

  ## plotting examples
  if(0){plot(out$betadraw)}
}
}

\keyword{datasets}
