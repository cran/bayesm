#include "bayesm.h"
 
// [[Rcpp::export]]
vec rdirichlet(vec const& alpha){
  
// Wayne Taylor 4/7/2015

// Purpose:
// draw from Dirichlet(alpha)

  int dim = alpha.size();
  vec y = zeros<vec>(dim);
  
  for(int i = 0; i<dim; i++) {    
      y[i] = rgamma(1,alpha[i])[0]; //rgamma returns a NumericVector, so adding [0] extracts the first element and treats it as type "double"
    }
  
  return(y/sum(y));
}
