#include "bayesm.h"
 
// [[Rcpp::export]]
double lndMvn(vec const& x, vec const& mu, mat const& rooti){

//Wayne Taylor 9/7/2014

// function to evaluate log of MV Normal density with  mean mu, var Sigma
// Sigma=t(root)%*%root   (root is upper tri cholesky root)
// Sigma^-1=rooti%*%t(rooti)   
// rooti is in the inverse of upper triangular chol root of sigma
//          note: this is the UL decomp of sigmai not LU!
//                Sigma=root'root   root=inv(rooti)

  vec z = vectorise(trans(rooti)*(x-mu));
  
  return((-(x.size()/2.0)*log(2*M_PI) -.5*(trans(z)*z) + sum(log(diagvec(rooti))))[0]);
}
