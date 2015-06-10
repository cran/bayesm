\name{margarine}
\alias{margarine}
\docType{data}
\title{Household Panel Data on Margarine Purchases}
\description{
  Panel data on purchases of margarine by 516 households. 
  Demographic variables are included.
}
\usage{data(margarine)}
\format{
  This is an R object that is a list of two data frames, list(choicePrice,demos)

 List of 2 \cr
 $ choicePrice:`data.frame':	4470 obs. of  12 variables:\cr
  \ldots $ hhid    : int [1:4470] 2100016 2100016 2100016 2100016 \cr
  \ldots $ choice  : num [1:4470] 1 1 1 1 1 4 1 1 4 1 \cr
  \ldots $ PPk_Stk : num [1:4470] 0.66 0.63 0.29 0.62 0.5 0.58 0.29  \cr
  \ldots $ PBB_Stk : num [1:4470] 0.67 0.67 0.5 0.61 0.58 0.45 0.51  \cr
  \ldots $ PFl_Stk : num [1:4470] 1.09 0.99 0.99 0.99 0.99 0.99 0.99 \cr
  \ldots $ PHse_Stk: num [1:4470] 0.57 0.57 0.57 0.57 0.45 0.45 0.29 \cr
  \ldots $ PGen_Stk: num [1:4470] 0.36 0.36 0.36 0.36 0.33 0.33 0.33 \cr
  \ldots $ PImp_Stk: num [1:4470] 0.93 1.03 0.69 0.75 0.72 0.72 0.72 \cr
  \ldots $ PSS_Tub : num [1:4470] 0.85 0.85 0.79 0.85 0.85 0.85 0.85 \cr
  \ldots $ PPk_Tub : num [1:4470] 1.09 1.09 1.09 1.09 1.07 1.07 1.07 \cr
  \ldots $ PFl_Tub : num [1:4470] 1.19 1.19 1.19 1.19 1.19 1.19 1.19 \cr
  \ldots $ PHse_Tub: num [1:4470] 0.33 0.37 0.59 0.59 0.59 0.59 0.59 \cr

  Pk is Parkay; BB is BlueBonnett, Fl is Fleischmanns, Hse is house,
  Gen is generic, Imp is Imperial, SS is Shed Spread.  _Stk indicates 
  stick, _Tub indicates Tub form.

 $ demos      :`data.frame':	516 obs. of  8 variables:\cr
  \ldots $ hhid     : num [1:516] 2100016 2100024 2100495 2100560 \cr
  \ldots $ Income   : num [1:516] 32.5 17.5 37.5 17.5 87.5 12.5 \cr
  \ldots $ Fs3_4    : int [1:516] 0 1 0 0 0 0 0 0 0 0 \cr
  \ldots $ Fs5      : int [1:516] 0 0 0 0 0 0 0 0 1 0 \cr
  \ldots $ Fam_Size : int [1:516] 2 3 2 1 1 2 2 2 5 2 \cr
  \ldots $ college  : int [1:516] 1 1 0 0 1 0 1 0 1 1 \cr
  \ldots $ whtcollar: int [1:516] 0 1 0 1 1 0 0 0 1 1 \cr
  \ldots $ retired  : int [1:516] 1 1 1 0 0 1 0 1 0 0 \cr

  Fs3_4 is dummy (family size 3-4). Fs5 is dummy for family size >= 5.
  college,whtcollar,retired are dummies reflecting these statuses.
}
\details{
  choice is a multinomial indicator of one of the 10 brands (in order listed under format). 
  All prices are in $.
}
\source{
  Allenby and Rossi (1991), "Quality Perceptions and Asymmetric Switching Between Brands," 
  \emph{Marketing Science} 10, 185-205.
}
\references{
  Chapter 5, \emph{Bayesian Statistics and Marketing} by Rossi et al.\cr
  \url{http://www.perossi.org/home/bsm-1}
}
\examples{
data(margarine)
cat(" Table of Choice Variable ",fill=TRUE)
print(table(margarine$choicePrice[,2]))
cat(" Means of Prices",fill=TRUE)
mat=apply(as.matrix(margarine$choicePrice[,3:12]),2,mean)
print(mat)
cat(" Quantiles of Demographic Variables",fill=TRUE)
mat=apply(as.matrix(margarine$demos[,2:8]),2,quantile)
print(mat)

##
## example of processing for use with rhierMnlRwMixture
##
if(0)
{
select= c(1:5,7)  ## select brands
chPr=as.matrix(margarine$choicePrice)
## make sure to log prices
chPr=cbind(chPr[,1],chPr[,2],log(chPr[,2+select]))
demos=as.matrix(margarine$demos[,c(1,2,5)])

## remove obs for other alts
chPr=chPr[chPr[,2] <= 7,]
chPr=chPr[chPr[,2] != 6,]

## recode choice
chPr[chPr[,2] == 7,2]=6

hhidl=levels(as.factor(chPr[,1]))
lgtdata=NULL
nlgt=length(hhidl)
p=length(select)  ## number of choice alts
ind=1
for (i in 1:nlgt) {
   nobs=sum(chPr[,1]==hhidl[i])
   if(nobs >=5) {
      data=chPr[chPr[,1]==hhidl[i],]
      y=data[,2]
      names(y)=NULL
      X=createX(p=p,na=1,Xa=data[,3:8],nd=NULL,Xd=NULL,INT=TRUE,base=1)
       lgtdata[[ind]]=list(y=y,X=X,hhid=hhidl[i]); ind=ind+1
   }
}
nlgt=length(lgtdata)
##
## now extract demos corresponding to hhs in lgtdata
##
Z=NULL
nlgt=length(lgtdata)
for(i in 1:nlgt){
   Z=rbind(Z,demos[demos[,1]==lgtdata[[i]]$hhid,2:3])
}
##
## take log of income and family size and demean
##
Z=log(Z)
Z[,1]=Z[,1]-mean(Z[,1])
Z[,2]=Z[,2]-mean(Z[,2])

keep=5
R=20000
mcmc1=list(keep=keep,R=R)
out=rhierMnlRwMixture(Data=list(p=p,lgtdata=lgtdata,Z=Z),Prior=list(ncomp=1),Mcmc=mcmc1)

summary(out$Deltadraw)
summary(out$nmix)

if(0){
## plotting examples
plot(out$nmix)
plot(out$Deltadraw)}
}

}
\keyword{datasets}
