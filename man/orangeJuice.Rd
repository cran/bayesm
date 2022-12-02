\name{orangeJuice}
\alias{orangeJuice}
\docType{data}

\title{Store-level Panel Data on Orange Juice Sales}

\description{
Weekly sales of refrigerated orange juice at 83 stores. Contains demographic information on those stores.
}

\usage{data(orangeJuice)}

\format{
The \code{orangeJuice} object is a list containing two data frames, \code{yx} and \code{storedemo}.
}

\details{
  In the \code{yx} data frame: 
  \tabular{ll}{
    \ldots\code{$store    } \tab store number \cr
    \ldots\code{$brand    } \tab brand indicator \cr
    \ldots\code{$week     } \tab week number \cr
    \ldots\code{$logmove  } \tab log of the number of units sold \cr
    \ldots\code{$constant } \tab a vector of 1s \cr
    \ldots\code{$price#   } \tab price of brand # \cr
    \ldots\code{$deal     } \tab in-store coupon activity \cr
    \ldots\code{$feature  } \tab feature advertisement \cr
    \ldots\code{$profit   } \tab profit 
  }
  The price variables correspond to the following brands:
  \tabular{ll}{
     1  \tab Tropicana Premium 64 oz \cr
     2  \tab Tropicana Premium 96 oz \cr
     3  \tab Florida's Natural 64 oz \cr
     4  \tab Tropicana 64 oz \cr
     5  \tab Minute Maid 64 oz \cr
     6  \tab Minute Maid 96 oz \cr
     7  \tab Citrus Hill 64 oz \cr
     8  \tab Tree Fresh 64 oz \cr
     9  \tab Florida Gold 64 oz \cr
     10 \tab Dominicks 64 oz \cr
     11 \tab Dominicks 128 oz 
  }
  In the \code{storedemo} data frame:
  \tabular{ll}{
    \ldots\code{$STORE    } \tab store number \cr
    \ldots\code{$AGE60    } \tab percentage of the population that is aged 60 or older \cr
    \ldots\code{$EDUC     } \tab percentage of the population that has a college degree \cr
    \ldots\code{$ETHNIC   } \tab percent of the population that is black or Hispanic \cr
    \ldots\code{$INCOME   } \tab median income \cr
    \ldots\code{$HHLARGE  } \tab percentage of households with 5 or more persons \cr
    \ldots\code{$WORKWOM  } \tab percentage of women with full-time jobs \cr
    \ldots\code{$HVAL150  } \tab percentage of households worth more than $150,000 \cr
    \ldots\code{$SSTRDIST } \tab distance to the nearest warehouse store \cr
    \ldots\code{$SSTRVOL  } \tab ratio of sales of this store to the nearest warehouse store \cr
    \ldots\code{$CPDIST5  } \tab average distance in miles to the nearest 5 supermarkets \cr
    \ldots\code{$CPWVOL5  } \tab ratio of sales of this store to the average of the nearest five stores 
  }
}

\source{Alan L. Montgomery (1997), "Creating Micro-Marketing Pricing Strategies Using Supermarket Scanner Data," \emph{Marketing Science} 16(4) 315--337.}

\references{Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
## load data
data(orangeJuice)

## print some quantiles of yx data  
cat("Quantiles of the Variables in yx data",fill=TRUE)
mat = apply(as.matrix(orangeJuice$yx), 2, quantile)
print(mat)

## print some quantiles of storedemo data
cat("Quantiles of the Variables in storedemo data",fill=TRUE)
mat = apply(as.matrix(orangeJuice$storedemo), 2, quantile)
print(mat)


## processing for use with rhierLinearModel
if(0) {
  
  ## select brand 1 for analysis
  brand1 = orangeJuice$yx[(orangeJuice$yx$brand==1),]
  
  store = sort(unique(brand1$store))
  nreg = length(store)
  nvar = 14
  
  regdata=NULL
  for (reg in 1:nreg) {
    y = brand1$logmove[brand1$store==store[reg]]
    iota = c(rep(1,length(y)))
    X = cbind(iota,log(brand1$price1[brand1$store==store[reg]]),
                   log(brand1$price2[brand1$store==store[reg]]),
                   log(brand1$price3[brand1$store==store[reg]]),
                   log(brand1$price4[brand1$store==store[reg]]),
                   log(brand1$price5[brand1$store==store[reg]]),
                   log(brand1$price6[brand1$store==store[reg]]),
                   log(brand1$price7[brand1$store==store[reg]]),
                   log(brand1$price8[brand1$store==store[reg]]),
                   log(brand1$price9[brand1$store==store[reg]]),
                   log(brand1$price10[brand1$store==store[reg]]),
                   log(brand1$price11[brand1$store==store[reg]]),
                   brand1$deal[brand1$store==store[reg]],
                   brand1$feat[brand1$store==store[reg]] )
    regdata[[reg]] = list(y=y, X=X)
    }
  
  ## storedemo is standardized to zero mean.
  Z = as.matrix(orangeJuice$storedemo[,2:12]) 
  dmean = apply(Z, 2, mean)
  for (s in 1:nreg) {Z[s,] = Z[s,] - dmean}
  iotaz = c(rep(1,nrow(Z)))
  Z = cbind(iotaz, Z)
  nz = ncol(Z)
  
  Data = list(regdata=regdata, Z=Z)
  Mcmc = list(R=R, keep=1)
  
  out = rhierLinearModel(Data=Data, Mcmc=Mcmc)
  
  summary(out$Deltadraw)
  summary(out$Vbetadraw)
  
  ## plotting examples
  if(0){ plot(out$betadraw) }
}
}

\keyword{datasets}


