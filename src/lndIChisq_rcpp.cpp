#include "bayesm.h"
 
// [[Rcpp::export]]
mat lndIChisq(double nu, double ssq, mat const& X) {

// Keunwoo Kim 07/24/2014

// Purpose: evaluate log-density of scaled Inverse Chi-sq density of random variable Z=nu*ssq/chisq(nu)
   
  return(-lgamma(nu/2)+(nu/2)*log((nu*ssq)/2)-((nu/2)+1)*log(X)-(nu*ssq)/(2*X));
}

