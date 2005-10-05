clusterMix=
function(zdraw,cutoff=.9,SILENT=FALSE){
#
#
# revision history:
#   written by p. rossi 9/05
#
# purpose: cluster observations based on draws of indicators of 
#   normal mixture components
#
# arguments:
#   zdraw is a R x nobs matrix of draws of indicators (typically output from rnmixGibbs)
#   the rth row of zdraw contains rth draw of indicators for each observations
#   each element of zdraw takes on up to p values for up to p groups. The maximum
#   number of groups is nobs.  Typically, however, the number of groups will be small
#   and equal to the number of components used in the normal mixture fit.
#
#   cutoff is a cutoff used in determining one clustering scheme it must be 
#   a number between .5 and 1.
#
# output:
#   two clustering schemes each with a vector of length nobs which gives the assignment
#   of each observation to a cluster
#
#   clustera (finds zdraw with similarity matrix closest to posterior mean of similarity)
#   clusterb (finds clustering scheme by assigning ones if posterior mean of similarity matrix
#             > cutoff and computing associated z )
#
# define needed functions
#
# ------------------------------------------------------------------------------------------   

ztoSim=function(z){
#
# function to convert indicator vector to Similarity matrix
# Sim is n x n matrix, Sim[i,j]=1 if pair(i,j) are in same group
# z is n x 1 vector of indicators (1,...,p)
#
# p.rossi 9/05
#
n=length(z)
zvec=c(rep(z,n))
zcomp=z%x%c(rep(1,n))
Sim=as.numeric((zvec==zcomp))
dim(Sim)=c(n,n)
return(Sim)
}
Simtoz=function(Sim){
#
# function to convert Similarity matrix to indicator vector
#  Sim is n x n matrix, Sim[i,j]=1 if pair(i,j) are in same group
#  z is vector of indicators from (1,...,p) of group memberships (dim n)
#
#
# p.rossi 9/05
n=ncol(Sim)
z=double(n)
i=1
groupn=1
while (i <= n){
  validind=z==0
  if(sum(Sim[validind,i]==1)>=1) {
     z[validind]=as.numeric(Sim[validind,i]==1)*groupn
     groupn=groupn+1
  }
  i=i+1
}
return(z)
} 
# ----------------------------------------------------------------------------------------
#
# check arguments
#
pandterm=function(message) { stop(message,call.=FALSE) }
if(missing(zdraw)) {pandterm("Requires zdraw argument -- R x n matrix of indicator draws")}
#
# check validity of zdraw rows -- must be integers in the range 1:nobs
#
nobs=ncol(zdraw)
R=nrow(zdraw)
if(sum(zdraw %in% (1:nobs)) < ncol(zdraw)*nrow(zdraw))
   {pandterm("Bad zdraw argument -- all elements must be integers in 1:nobs")}
cat("Table of zdraw values pooled over all rows",fill=TRUE)
print(table(zdraw))
#
# check validity of cuttoff
if(cutoff > 1 || cutoff < .5) {pandterm(paste("cutoff invalid, = ",cutoff))}

#
# compute posterior mean of Similarity matrix
#
#
if(!SILENT){
   cat("Computing Posterior Expectation of Similarity Matrix",fill=TRUE)
   cat("processing draws ...",fill=TRUE); fsh()
}
Pmean=matrix(0,nrow=nobs,ncol=nobs)
R=nrow(zdraw)
for (r in 1:R) {
   Pmean=Pmean+ztoSim(out$zdraw[r,])
   if(!SILENT) {if(r%%100 == 0) {cat("  ",r,fill=TRUE); fsh()}}
}
Pmean=Pmean/R

#
# now find index for draw which minimizes discrepancy between
# post exp of similarity and sim implied by that z
if(!SILENT){
  cat(" ",fill=TRUE)
  cat("Look for zdraw which minimizes loss",fill=TRUE)
  cat("processing draws ...",fill=TRUE); fsh()
}
loss=double(R)
for (r in 1:R){
  loss[r]=sum(abs(Pmean-ztoSim(out$zdraw[r,]))) 
  if(!SILENT) {if(r%%100 == 0) {cat("  ",r,fill=TRUE);fsh()}}
}
index=which(loss==min(loss))
clustera=zdraw[index[1],]
#
# now due clustering by assigning Similarity to any (i,j) pair for which
# Pmean > cutoff
Sim=matrix(as.numeric(Pmean >= cutoff),ncol=nobs)
clusterb=Simtoz(Sim)
return(list(clustera=clustera,clusterb=clusterb))
}
   
      



