#include "bayesm.h"
 
// [[Rcpp::export]]
NumericVector rtrun(NumericVector const& mu, NumericVector const& sigma, 
                         NumericVector const& a, NumericVector const& b){
                           
// Wayne Taylor 9/7/2014

// function to draw from univariate truncated norm
// a is vector of lower bounds for truncation
// b is vector of upper bounds for truncation

  NumericVector FA = pnorm((a-mu)/sigma);
  NumericVector FB = pnorm((b-mu)/sigma);
  
  return(mu+sigma*qnorm(runif(mu.size())*(FB-FA)+FA));
}
