mixDenBi=
function(i,j,xi,xj,pvec,comps) 
{
# Revision History:
#   P. Rossi 6/05
#
# purpose: compute marg bivariate density implied by mixture of multivariate normals specified
#			by pvec,comps
#
# arguments:
#     i,j:  index of two variables
#     xi specifies a grid of points for var i
#     xj specifies a grid of points for var j
#     pvec: prior probabilities of normal components
#     comps: list, each member is a list comp with ith normal component ~ N(comp[[1]],Sigma), 
#            Sigma = t(R)%*%R, R^{-1} = comp[[2]]
# Output:
#     matrix with values of density on grid
#
# ---------------------------------------------------------------------------------------------
# define function needed
#
bivcomps=function(i,j,comps)
{
# purpose: obtain marginal means and standard deviations from list of normal components
# arguments:
#     i,j: index of elements for bivariate marginal
#     comps: list, each member is a list comp with ith normal component ~N(comp[[1]],Sigma), 
#            Sigma = t(R)%*%R, R^{-1} = comp[[2]]
# returns:
#  a list with relevant mean vectors and rooti for each compenent
#  [[2]]=$sigma a matrix whose ith row is the standard deviations for the ith component
#
result=NULL
nc = length(comps)
dim = length(comps[[1]][[1]])
ind=matrix(c(i,j,i,j,i,i,j,j),ncol=2)
for(comp in 1:nc) {
   mu = comps[[comp]][[1]][c(i,j)]
   root= backsolve(comps[[comp]][[2]],diag(dim))
   Sigma=crossprod(root)
   sigma=matrix(Sigma[ind],ncol=2)
   rooti=backsolve(chol(sigma),diag(2))
   result[[comp]]=list(mu=mu,rooti=rooti)
}
return(result)
}
normden=
function(x,mu,rooti)
{
#
# function to evaluate MV NOrmal density with  mean mu, var Sigma
# and with sigma^-1=rooti%*%t(rooti)   
# rooti is in the inverse of upper triangular chol root of sigma
#          note: this is the UL decomp of sigmai not LU!
#                Sigma=root'root   root=inv(rooti)
#
z=as.vector(t(rooti)%*%(x-mu))
exp(-.5*length(x)*log(2*pi)-.5*(z%*%z) + sum(log(diag(rooti))))
}
# ----------------------------------------------------------------------------------------------
nc = length(comps)
marmoms=bivcomps(i,j,comps)
den = matrix(0.0,nrow=length(xi),ncol=length(xj))
for(indi in 1:length(xi)) {
   for(indj in 1:length(xj)) {
      for(comp in 1:nc) {
          den[indi,indj] = den[indi,indj] + normden(c(xi[indi],xj[indj]),marmoms[[comp]]$mu,
                                marmoms[[comp]]$rooti)*pvec[comp]
      }
    } 
}
return(den)
}
