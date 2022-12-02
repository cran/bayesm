\name{margarine}
\alias{margarine}
\docType{data}

\title{Household Panel Data on Margarine Purchases}

\description{
  Panel data on purchases of margarine by 516 households. Demographic variables are included.
}

\usage{data(margarine)}

\format{
The \code{detailing} object is a list containing two data frames, \code{choicePrice} and \code{demos}.
}

\details{
In the \code{choicePrice} data frame:
  \tabular{ll}{
    \ldots\code{$hhid  } \tab household ID \cr
    \ldots\code{$choice} \tab multinomial indicator of one of the 10 products 
  }
The products are indicated by brand and type. 

Brands:
  \tabular{ll}{
    \ldots\code{$Pk  } \tab Parkay \cr
    \ldots\code{$BB  } \tab BlueBonnett \cr
    \ldots\code{$Fl  } \tab Fleischmanns \cr
    \ldots\code{$Hse } \tab house \cr
    \ldots\code{$Gen } \tab generic \cr
    \ldots\code{$Imp } \tab Imperial \cr
    \ldots\code{$SS  } \tab Shed Spread 
  }
Product type:
  \tabular{ll}{
    \ldots\code{$_Stk} \tab stick \cr
    \ldots\code{$_Tub} \tab tub 
  }  
In the \code{demos} data frame:
  \tabular{ll}{
    \ldots\code{$Fs3_4    } \tab dummy for family size 3-4 \cr
    \ldots\code{$Fs5      } \tab dummy for family size >= 5 \cr
    \ldots\code{$college  } \tab dummy for education status \cr
    \ldots\code{$whtcollar} \tab dummy for job status \cr
    \ldots\code{$retired  } \tab dummy for retirement status 
  }
All prices are in U.S. dollars.
}

\source{Allenby, Greg and Peter Rossi (1991), "Quality Perceptions and Asymmetric Switching Between Brands," \emph{Marketing Science} 10, 185--205.
}

\references{
Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi, Allenby, and McCulloch.}

\examples{
data(margarine)
cat(" Table of Choice Variable ", fill=TRUE)
print(table(margarine$choicePrice[,2]))

cat(" Means of Prices", fill=TRUE)
mat=apply(as.matrix(margarine$choicePrice[,3:12]), 2, mean)
print(mat)

cat(" Quantiles of Demographic Variables", fill=TRUE)
mat=apply(as.matrix(margarine$demos[,2:8]), 2, quantile)
print(mat)


## example of processing for use with 'rhierMnlRwMixture'
if(0) {
  select = c(1:5,7)  ## select brands
  chPr = as.matrix(margarine$choicePrice)
  
  ## make sure to log prices
  chPr = cbind(chPr[,1], chPr[,2], log(chPr[,2+select]))
  demos = as.matrix(margarine$demos[,c(1,2,5)])
  
  ## remove obs for other alts
  chPr = chPr[chPr[,2] <= 7,]
  chPr = chPr[chPr[,2] != 6,]
  
  ## recode choice
  chPr[chPr[,2] == 7,2] = 6
  
  hhidl = levels(as.factor(chPr[,1]))
  lgtdata = NULL
  nlgt = length(hhidl)
  p = length(select)  ## number of choice alts
  
  ind = 1
  for (i in 1:nlgt) {
    nobs = sum(chPr[,1]==hhidl[i])
    if(nobs >=5) {
      data = chPr[chPr[,1]==hhidl[i],]
      y = data[,2]
      names(y) = NULL
      X = createX(p=p, na=1, Xa=data[,3:8], nd=NULL, Xd=NULL, INT=TRUE, base=1)
      lgtdata[[ind]] = list(y=y, X=X, hhid=hhidl[i])
      ind = ind+1
    }
  }
  nlgt = length(lgtdata)
  
  ## now extract demos corresponding to hhs in lgtdata
  Z = NULL
  nlgt = length(lgtdata)
  for(i in 1:nlgt){
     Z = rbind(Z, demos[demos[,1]==lgtdata[[i]]$hhid, 2:3])
  }
  
  ## take log of income and family size and demean
  Z = log(Z)
  Z[,1] = Z[,1] - mean(Z[,1])
  Z[,2] = Z[,2] - mean(Z[,2])
  
  keep = 5
  R = 20000
  mcmc1 = list(keep=keep, R=R)
  
  out = rhierMnlRwMixture(Data=list(p=p,lgtdata=lgtdata, Z=Z),
                          Prior=list(ncomp=1), Mcmc=mcmc1)
  
  summary(out$Deltadraw)
  summary(out$nmix)
  
  ## plotting examples
  if(0){
    plot(out$nmix)
    plot(out$Deltadraw)
  }
}
}

\keyword{datasets}
