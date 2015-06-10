\name{detailing}
\alias{detailing}
\docType{data}
\title{ Physician Detailing Data from Manchanda et al (2004)}
\description{
  Monthly data on detailing (sales calls) on 1000 physicians. 23 mos of data
  for each physician. Includes physician covariates. Dependent variable (\code{scripts}) is the
  number of new prescriptions ordered by the physician for the drug detailed.
}
\usage{data(detailing)}
\format{
 This R object is a list of two data frames, list(counts,demo).
 
 List of 2:

 $ counts:`data.frame':	23000 obs. of  4 variables:\cr
  \ldots$ id        : int [1:23000] 1 1 1 1 1 1 1 1 1 1 \cr
  \ldots$ scripts       : int [1:23000] 3 12 3 6 5 2 5 1 5 3 \cr
  \ldots$ detailing     : int [1:23000] 1 1 1 2 1 0 2 2 1 1 \cr
  \ldots$ lagged_scripts: int [1:23000] 4 3 12 3 6 5 2 5 1 5 

 $ demo  :`data.frame':	1000 obs. of  4 variables:\cr
  \ldots$ id          : int [1:1000] 1 2 3 4 5 6 7 8 9 10 \cr
  \ldots$ generalphys : int [1:1000] 1 0 1 1 0 1 1 1 1 1 \cr
  \ldots$ specialist: int [1:1000] 0 1 0 0 1 0 0 0 0 0  \cr
  \ldots$ mean_samples: num [1:1000] 0.722 0.491 0.339 3.196 0.348 
}
\details{
  generalphys is dummy for if doctor is a "general practitioner," specialist is dummy for
  if the physician is a specialist in the theraputic class for which the drug is 
  intended, mean_samples is the mean number of free drug samples given the doctor
  over the sample.
}
\source{
  Manchanda, P., P. K. Chintagunta and P. E. Rossi (2004), "Response Modeling with Non-Random
  Marketing Mix Variables," \emph{Journal of Marketing Research} 41, 467-478.
}
\examples{
data(detailing)
cat(" table of Counts Dep Var", fill=TRUE)
print(table(detailing$counts[,2]))
cat(" means of Demographic Variables",fill=TRUE)
mat=apply(as.matrix(detailing$demo[,2:4]),2,mean)
print(mat)

##
## example of processing for use with rhierNegbinRw
##
if(0)
{
data(detailing)
counts = detailing$counts
Z = detailing$demo

# Construct the Z matrix
Z[,1] = 1
Z[,2]=Z[,2]-mean(Z[,2])
Z[,3]=Z[,3]-mean(Z[,3])
Z[,4]=Z[,4]-mean(Z[,4])
Z=as.matrix(Z)
id=levels(factor(counts$id))
nreg=length(id)
nobs = nrow(counts$id)

regdata=NULL
for (i in 1:nreg) {
    X = counts[counts[,1] == id[i],c(3:4)]
    X = cbind(rep(1,nrow(X)),X)
    y = counts[counts[,1] == id[i],2]
    X = as.matrix(X)
    regdata[[i]]=list(X=X, y=y)
}
nvar=ncol(X)            # Number of X variables
nz=ncol(Z)              # Number of Z variables
rm(detailing,counts)              
cat("Finished Reading data",fill=TRUE)
fsh()

Data = list(regdata=regdata, Z=Z)
deltabar = matrix(rep(0,nvar*nz),nrow=nz)
Vdelta = 0.01 * diag(nz)
nu = nvar+3
V = 0.01*diag(nvar)
a = 0.5
b = 0.1
Prior = list(deltabar=deltabar, Vdelta=Vdelta, nu=nu, V=V, a=a, b=b)

R = 10000
keep =1
s_beta=2.93/sqrt(nvar)
s_alpha=2.93
c=2
Mcmc = list(R=R, keep = keep, s_beta=s_beta, s_alpha=s_alpha, c=c)
out = rhierNegbinRw(Data, Prior, Mcmc)

# Unit level mean beta parameters
Mbeta = matrix(rep(0,nreg*nvar),nrow=nreg)
ndraws = length(out$alphadraw)
for (i in 1:nreg) { Mbeta[i,] = rowSums(out$Betadraw[i, , ])/ndraws }

cat(" Deltadraws ",fill=TRUE)
summary(out$Deltadraw)
cat(" Vbetadraws ",fill=TRUE)
summary(out$Vbetadraw)
cat(" alphadraws ",fill=TRUE)
summary(out$alphadraw)

if(0){
## plotting examples
plot(out$betadraw)
plot(out$alphadraw)
plot(out$Deltadraw)
}
}

}
\keyword{datasets}
