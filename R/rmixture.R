rmixture=
function(n,p,comps)
{
#
# R. McCulloch 12/04
# revision history:
#   commented by rossi 3/05
#
# purpose: iid draws from mixture of multivariate normals
# arguments:
#     n: number of draws
#     p: prior probabilities of normal components
#     comps: list, each member is a list comp with ith normal component
#                     ~N(comp[[1]],Sigma), Sigma = t(R)%*%R, R^{-1} = comp[[2]]
# output:
#  list of x (n by length(comp[[1]]) matrix of draws) and z latent indicators of
#  component
#
#----------------------------------------------------------------------------------
# define function needed
#
rcomp=function(comp) {
# purpose: draw multivariate normal with mean and variance given by comp 
# arguments:
#     comp is a list of length 2,
#     comp[[1]] is the mean and comp[[2]] is R^{-1} = comp[[2]], Sigma = t(R)%*%R
invUT = function(ut) {
backsolve(ut,diag(rep(1,nrow(ut))))
}
as.vector(comp[[1]] + t(invUT(comp[[2]]))%*%rnorm(length(comp[[1]])))
}
#----------------------------------------------------------------------------------
#
z = sample(1:length(p), n, replace = TRUE, prob = p)
return(list(x = t(sapply(comps[z],rcomp)),z=z))
}
