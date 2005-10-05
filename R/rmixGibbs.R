rmixGibbs=
function(y,Bbar,A,nu,V,a,p,z,comps)
{
#
# Revision History: 
#   R. McCulloch 11/04
#   P. Rossi 3/05 put in backsolve and improved documentation
#
# purpose: do gibbs sampling inference for a mixture of multivariate normals
#
# arguments:
#     y: data, rows are observations, assumed to be iid draws from normal mixture
#     Bbar,A,nu,V: common prior for mean and variance of each normal component
#
#     note: Bbar should be a matrix. usually with only one row
#
#        beta ~ N(betabar,Sigma (x) A^-1)
#                   betabar=vec(Bbar)
#          Sigma ~ IW(nu,V) or Sigma^-1 ~ W(nu, V^-1)
#              note: if you want Sigma ~ A, use nu big and rwishart(nu,nu(A)^{-1})$IW
#     a: Dirichlet parameters for prior on p
#     p: prior probabilities of normal components
#     z: components indentities for each observation 
#        (vector of intergers each in {1,2,...number of components})
#     comps: list, each member is a list comp with ith normal component
#                        ~N(comp[[1]],Sigma), Sigma = t(R)%*%R, R^{-1} = comp[[2]]
# Output:
#  list with elements [[1]=$p, [[2]]=$z, and [[3]]=$comps, with the updated values
#
#------------------------------------------------------------------------------------
#  define functions needed
#
rcompsC = function(x,p,comps) {
# purpose:
#  draws class membership of rows of x, given x rows are iid draws from 
#  mixture of multivariate normals
# arguments:
#     x: observations (number of observations x dimension)
#     p: prior probabilities of mixture components
#     comps: list, each member is a list with mean and R^{-1}, Sigma = t(R)%*%R 
dim = ncol(x)
nob = nrow(x)
nc = length(comps)
mumat = matrix(0.0,dim,nc)
rivmat = matrix(0.0,dim*(dim+1)/2,nc)
for(i in 1:nc) {
   mumat[,i] = comps[[i]][[1]]
   rivmat[,i] = uttovC(comps[[i]][[2]])
}
xx=t(x)
.C('crcomps',as.double(xx),as.double(mumat),as.double(rivmat),as.double(p),
    as.integer(dim),as.integer(nc),as.integer(nob),res=integer(nob))$res
}

uttovC = function(rooti) {
# returns vector of square upper triangular matrix rooti, goes down columns dropping the zeros
dim = nrow(rooti)
n = dim*(dim+1)/2
.C('cuttov',as.double(rooti),res = double(n),as.integer(dim))$res
}
#-----------------------------------------------------------------------------------------
nmix = length(a)
#draw comps
for(i in 1:nmix) {
   nobincomp = sum(z==i)         # get number of observations "in" component i
   if(nobincomp>0) {             # if more than one obs in this component, draw from posterior
      yi=y[z==i,]
      dim(yi)=c(nobincomp,ncol(y))
          #  worry about case where y has only one col (univ mixtures) or only one row
          #  then yi gets converted to a vector
      temp = rmultireg(yi,matrix(rep(1,nobincomp),ncol=1),Bbar,A,nu,V)
      comps[[i]] = list(mu = as.vector(temp$B),
                        rooti=backsolve(chol(temp$Sigma),diag(rep(1,nrow(temp$Sigma)))))
   } 
   else { # else draw from the prior
      rw=rwishart(nu,chol2inv(chol(V)))
      comps[[i]] = list(mu = as.vector(t(Bbar) + (rw$CI %*% rnorm(length(Bbar)))/sqrt(A[1,1])),
                        rooti=backsolve(chol(rw$IW),diag(rep(1,nrow(V)))))
   }
}
#draw z
z=rcompsC(y,p,comps)
#draw p
for(i in 1:length(a)) a[i] = a[i] + sum(z==i)
p = rdirichlet(a)
return(list(p=p,z=z,comps=comps))
}
