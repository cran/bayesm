rmvst=
function(nu,mu,root){
#
# function to draw from MV s-t  with nu df, mean mu, Sigma=t(root)%*%root
#      root is upper triangular cholesky root
nvec=t(root)%*%rnorm(length(mu))
return(nvec/sqrt(rchisq(1,nu)/nu) + mu)
}
