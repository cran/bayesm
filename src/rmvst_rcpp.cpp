#include "bayesm.h"
 
// [[Rcpp::export]]
vec rmvst(int nu, vec const& mu, mat const& root){
  
// Wayne Taylor 9/7/2014

// function to draw from MV s-t  with nu df, mean mu, Sigma=t(root)%*%root
//  root is upper triangular cholesky root

  vec rnormd = rnorm(mu.size());
  vec nvec = trans(root)*rnormd;
  
  return(nvec/sqrt(rchisq(1,nu)[0]/nu) + mu); //rchisq returns a vectorized object, so using [0] allows for the conversion to double
}
