\name{bank}
\alias{bank}
\docType{data}
\title{ Bank Card Conjoint Data of Allenby and Ginter (1995)}
\description{
  Data from a conjoint experiment in which two partial profiles of 
  credit cards were presented to 946 respondents. The variable 
  bank$choiceAtt$choice indicates which profile was chosen.  The
  profiles are coded as the difference in attribute levels. Thus,
  a "-1" means the profile coded as a choice of "0" has the attribute.
  A value of 0 means that the attribute was not present in the 
  comparison.
  
  data on age,income and gender (female=1) are also recorded in 
  bank$demo
}
\usage{data(bank)}
\format{
  This R object is a list of two data frames, list(choiceAtt,demo).

  List of 2 

 $ choiceAtt:`data.frame':	14799 obs. of  16 variables:\cr
  \ldots$ id           : int [1:14799] 1 1 1 1 1 1 1 1 1 1 \cr
  \ldots$ choice       : int [1:14799] 1 1 1 1 1 1 1 1 0 1 \cr
  \ldots$ Med_FInt     : int [1:14799] 1 1 1 0 0 0 0 0 0 0 \cr
  \ldots$ Low_FInt     : int [1:14799] 0 0 0 0 0 0 0 0 0 0 \cr
  \ldots$ Med_VInt     : int [1:14799] 0 0 0 0 0 0 0 0 0 0 \cr
  \ldots$ Rewrd_2      : int [1:14799] -1 1 0 0 0 0 0 1 -1 0 \cr
  \ldots$ Rewrd_3      : int [1:14799] 0 -1 1 0 0 0 0 0 1 -1 \cr
  \ldots$ Rewrd_4      : int [1:14799] 0 0 -1 0 0 0 0 0 0 1 \cr
  \ldots$ Med_Fee      : int [1:14799] 0 0 0 1 1 -1 -1 0 0 0 \cr
  \ldots$ Low_Fee      : int [1:14799] 0 0 0 0 0 1 1 0 0 0 \cr
  \ldots$ Bank_B       : int [1:14799] 0 0 0 -1 1 -1 1 0 0 0 \cr
  \ldots$ Out_State    : int [1:14799] 0 0 0 0 -1 0 -1 0 0 0 \cr
  \ldots$ Med_Rebate   : int [1:14799] 0 0 0 0 0 0 0 0 0 0 \cr
  \ldots$ High_Rebate  : int [1:14799] 0 0 0 0 0 0 0 0 0 0 \cr
  \ldots$ High_CredLine: int [1:14799] 0 0 0 0 0 0 0 -1 -1 -1 \cr
  \ldots$ Long_Grace   : int [1:14799] 0 0 0 0 0 0 0 0 0 0 

 $ demo     :`data.frame':	946 obs. of  4 variables:\cr
  \ldots$ id    : int [1:946] 1 2 3 4 6 7 8 9 10 11 \cr
  \ldots$ age   : int [1:946] 60 40 75 40 30 30 50 50 50 40 \cr
  \ldots$ income: int [1:946] 20 40 30 40 30 60 50 100 50 40 \cr
  \ldots$ gender: int [1:946] 1 1 0 0 0 0 1 0 0 0 \cr
}
\details{
  Each respondent was presented with between 13 and 17 paired comparisons. Thus, this
  dataset has a panel structure.
}
\source{
  Allenby and Ginter (1995), "Using Extremes to Design Products and Segment
  Markets," \emph{JMR}, 392-403.
}
\references{ Appendix A, \emph{Bayesian Statistics and Marketing}
  by Rossi, Allenby and McCulloch. \cr
  \url{http://www.perossi.org/home/bsm-1}
}
\examples{
data(bank)
cat(" table of Binary Dep Var", fill=TRUE)
print(table(bank$choiceAtt[,2]))
cat(" table of Attribute Variables",fill=TRUE)
mat=apply(as.matrix(bank$choiceAtt[,3:16]),2,table)
print(mat)
cat(" means of Demographic Variables",fill=TRUE)
mat=apply(as.matrix(bank$demo[,2:3]),2,mean)
print(mat)

## example of processing for use with rhierBinLogit
##
if(0)
{
choiceAtt=bank$choiceAtt
Z=bank$demo

## center demo data so that mean of random-effects
## distribution can be interpreted as the average respondent

Z[,1]=rep(1,nrow(Z))
Z[,2]=Z[,2]-mean(Z[,2])
Z[,3]=Z[,3]-mean(Z[,3])
Z[,4]=Z[,4]-mean(Z[,4])
Z=as.matrix(Z)

hh=levels(factor(choiceAtt$id))
nhh=length(hh)
lgtdata=NULL
for (i in 1:nhh) {
	y=choiceAtt[choiceAtt[,1]==hh[i],2]
	nobs=length(y)
	X=as.matrix(choiceAtt[choiceAtt[,1]==hh[i],c(3:16)])
	lgtdata[[i]]=list(y=y,X=X)
		}

cat("Finished Reading data",fill=TRUE)
fsh()

Data=list(lgtdata=lgtdata,Z=Z)
Mcmc=list(R=10000,sbeta=0.2,keep=20)
set.seed(66)
out=rhierBinLogit(Data=Data,Mcmc=Mcmc)

begin=5000/20
end=10000/20

summary(out$Deltadraw,burnin=begin)
summary(out$Vbetadraw,burnin=begin)

if(0){
## plotting examples

## plot grand means of random effects distribution (first row of Delta)
index=4*c(0:13)+1
matplot(out$Deltadraw[,index],type="l",xlab="Iterations/20",ylab="",
main="Average Respondent Part-Worths")

## plot hierarchical coefs
plot(out$betadraw)

## plot log-likelihood
plot(out$llike,type="l",xlab="Iterations/20",ylab="",main="Log Likelihood")

}
}
}
\keyword{datasets}
