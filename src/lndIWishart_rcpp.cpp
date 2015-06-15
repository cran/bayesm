#include "bayesm.h"
 
// [[Rcpp::export]]
double lndIWishart(double nu, mat const& V, mat const& IW){

// Keunwoo Kim 07/24/2014

// Purpose: evaluate log-density of inverted Wishart with normalizing constant

// Arguments: 
//        nu is d. f. parm
//        V is location matrix
//        IW is the value at which the density should be evaluated

// Note: in this parameterization, E[IW]=V/(nu-k-1)

  int k = V.n_cols;
  mat Uiw = chol(IW);
  mat Uiwi = solve(trimatu(Uiw), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  mat IWi = Uiwi*trans(Uiwi);
  mat cholV = chol(V);
  double lndetVd2 = sum(log(cholV.diag()));
  double lndetIWd2 = sum(log(Uiw.diag()));
  
  // first evaluate constant
  double cnst = ((nu*k)/2)*log(2.0)+((k*(k-1))/4.0)*log(M_PI); // (k*(k-1))/4 is recognized as integer. "4.0" allows it to be recognized as a double.
  vec seq_1_k = cumsum(ones<vec>(k)); // build c(1:k) through cumsum function
  vec arg = (nu+1-seq_1_k)/2.0;
  
  // lgamma cannot receive arma::vec input. Compute cnst+sum(lgamma(arg)).
  for (int i=0; i<k; i++){
    cnst = cnst+lgamma(arg[i]);
  }
  
  return (-cnst+nu*lndetVd2-(nu+k+1)*lndetIWd2-.5*sum((V*IWi.diag())));
}
