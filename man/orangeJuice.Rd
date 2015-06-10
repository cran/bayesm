\name{orangeJuice}
\alias{orangeJuice}
\docType{data}
\title{Store-level Panel Data on Orange Juice Sales}
\description{
  yx, weekly sales of refrigerated orange juice at 83 stores. \cr
  storedemo, contains demographic information on those stores. \cr
}
\usage{data(orangeJuice)}
\format{
 This R object is a list of two data frames, list(yx,storedemo).\cr

 List of 2 \cr
 $ yx       :'data.frame':	106139 obs. of  19 variables:\cr
  \ldots $ store   : int [1:106139] 2 2 2 2 2 2 2 2 2 2 \cr
  \ldots $ brand   : int [1:106139] 1 1 1 1 1 1 1 1 1 1 \cr
  \ldots $ week    : int [1:106139] 40 46 47 48 50 51 52 53 54 57 \cr
  \ldots $ logmove : num [1:106139] 9.02 8.72 8.25 8.99 9.09 \cr
  \ldots $ constant: int [1:106139] 1 1 1 1 1 1 1 1 1 1 \cr
  \ldots $ price1  : num [1:106139] 0.0605 0.0605 0.0605 0.0605 0.0605 \cr
  \ldots $ price2  : num [1:106139] 0.0605 0.0603 0.0603 0.0603 0.0603 \cr
  \ldots $ price3  : num [1:106139] 0.0420 0.0452 0.0452 0.0498 0.0436 \cr
  \ldots $ price4  : num [1:106139] 0.0295 0.0467 0.0467 0.0373 0.0311 \cr
  \ldots $ price5  : num [1:106139] 0.0495 0.0495 0.0373 0.0495 0.0495 \cr
  \ldots $ price6  : num [1:106139] 0.0530 0.0478 0.0530 0.0530 0.0530 \cr
  \ldots $ price7  : num [1:106139] 0.0389 0.0458 0.0458 0.0458 0.0466 \cr
  \ldots $ price8  : num [1:106139] 0.0414 0.0280 0.0414 0.0414 0.0414 \cr
  \ldots $ price9  : num [1:106139] 0.0289 0.0430 0.0481 0.0423 0.0423 \cr
  \ldots $ price10 : num [1:106139] 0.0248 0.0420 0.0327 0.0327 0.0327 \cr
  \ldots $ price11 : num [1:106139] 0.0390 0.0390 0.0390 0.0390 0.0382 \cr
  \ldots $ deal    : int [1:106139] 1 0 0 0 0 0 1 1 1 1 \cr
  \ldots $ feat    : num [1:106139] 0 0 0 0 0 0 0 0 0 0 \cr
  \ldots $ profit  : num [1:106139] 38.0 30.1 30.0 29.9 29.9 \cr

     1 Tropicana Premium 64 oz;   2 Tropicana Premium 96 oz;  3 Florida's Natural 64 oz; \cr   
     4 Tropicana 64 oz;           5 Minute Maid 64 oz;        6 Minute Maid 96 oz; \cr
     7 Citrus Hill 64 oz;         8 Tree Fresh 64 oz;         9 Florida Gold 64 oz; \cr       
     10 Dominicks 64 oz;          11 Dominicks 128 oz.  \cr

 $ storedemo:'data.frame':	83 obs. of  12 variables:\cr
  \ldots $ STORE   : int [1:83] 2 5 8 9 12 14 18 21 28 32 \cr
  \ldots $ AGE60   : num [1:83] 0.233 0.117 0.252 0.269 0.178 \cr
  \ldots $ EDUC    : num [1:83] 0.2489 0.3212 0.0952 0.2222 0.2534 \cr
  \ldots $ ETHNIC  : num [1:83] 0.1143 0.0539 0.0352 0.0326 0.3807 \cr
  \ldots $ INCOME  : num [1:83] 10.6 10.9 10.6 10.8 10.0 \cr
  \ldots $ HHLARGE : num [1:83] 0.1040 0.1031 0.1317 0.0968 0.0572 \cr
  \ldots $ WORKWOM : num [1:83] 0.304 0.411 0.283 0.359 0.391 \cr
  \ldots $ HVAL150 : num [1:83] 0.4639 0.5359 0.0542 0.5057 0.3866 \cr
  \ldots $ SSTRDIST: num [1:83] 2.11 3.80 2.64 1.10 9.20 \cr
  \ldots $ SSTRVOL : num [1:83] 1.143 0.682 1.500 0.667 1.111 \cr
  \ldots $ CPDIST5 : num [1:83] 1.93 1.60 2.91 1.82 0.84 \cr
  \ldots $ CPWVOL5 : num [1:83] 0.377 0.736 0.641 0.441 0.106 \cr
}
\details{
  \describe{
    \item{\code{store}}{store number}
    \item{\code{brand}}{brand indicator}
    \item{\code{week}}{week number}
    \item{\code{logmove}}{log of the number of units sold}
    \item{\code{constant}}{a vector of 1}
    \item{\code{price1}}{price of brand 1}
    \item{\code{deal}}{in-store coupon activity}
    \item{\code{feature}}{feature advertisement}
    \item{\code{STORE}}{store number}
    \item{\code{AGE60}}{percentage of the population that is aged 60 or older}
    \item{\code{EDUC}}{percentage of the population that has a college degree}
    \item{\code{ETHNIC}}{percent of the population that is black or Hispanic}
    \item{\code{INCOME}}{median income}
    \item{\code{HHLARGE}}{percentage of households with 5 or more persons}
    \item{\code{WORKWOM}}{percentage of women with full-time jobs}
    \item{\code{HVAL150}}{percentage of households worth more than $150,000}
    \item{\code{SSTRDIST}}{distance to the nearest warehouse store}
    \item{\code{SSTRVOL}}{ratio of sales of this store to the nearest warehouse store}
    \item{\code{CPDIST5}}{average distance in miles to the nearest 5 supermarkets}
    \item{\code{CPWVOL5}}{ratio of sales of this store to the average of the nearest five stores}
    }
}
\source{
   Alan L. Montgomery (1997), "Creating Micro-Marketing Pricing Strategies Using Supermarket Scanner Data,"
   \emph{Marketing Science} 16(4) 315-337. 
}
\references{
    Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi et al.\cr
    \url{http://www.perossi.org/home/bsm-1}
}
\examples{

## Example 
## load data
data(orangeJuice)

## print some quantiles of yx data  
cat("Quantiles of the Variables in yx data",fill=TRUE)
mat=apply(as.matrix(orangeJuice$yx),2,quantile)
print(mat)

## print some quantiles of storedemo data
cat("Quantiles of the Variables in storedemo data",fill=TRUE)
mat=apply(as.matrix(orangeJuice$storedemo),2,quantile)
print(mat)


## Example 2 processing for use with rhierLinearModel
##
##
if(0)
{

## select brand 1 for analysis
brand1=orangeJuice$yx[(orangeJuice$yx$brand==1),]

store = sort(unique(brand1$store))
nreg = length(store)
nvar=14

regdata=NULL
for (reg in 1:nreg) {
        y=brand1$logmove[brand1$store==store[reg]]
        iota=c(rep(1,length(y)))
        X=cbind(iota,log(brand1$price1[brand1$store==store[reg]]),
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
                     brand1$feat[brand1$store==store[reg]])
        regdata[[reg]]=list(y=y,X=X)
      }

## storedemo is standardized to zero mean.

Z=as.matrix(orangeJuice$storedemo[,2:12]) 
dmean=apply(Z,2,mean)
for (s in 1:nreg){
        Z[s,]=Z[s,]-dmean
}
iotaz=c(rep(1,nrow(Z)))
Z=cbind(iotaz,Z)
nz=ncol(Z)


Data=list(regdata=regdata,Z=Z)
Mcmc=list(R=R,keep=1)

out=rhierLinearModel(Data=Data,Mcmc=Mcmc)

summary(out$Deltadraw)
summary(out$Vbetadraw)

if(0){
## plotting examples
plot(out$betadraw)
}
}

}
\keyword{datasets}


