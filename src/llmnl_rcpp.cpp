#include "bayesm.h"
 
//[[Rcpp::export]]
double llmnl(vec const& beta, vec const& y, mat const& X){
  
// Wayne Taylor 9/7/2014

// Evaluates log-likelihood for the multinomial logit model

  int n = y.size();
  int j = X.n_rows/n;
  mat Xbeta = X*beta;
      
  vec xby = zeros<vec>(n);
  vec denom = zeros<vec>(n);
  
  for(int i = 0; i<n;i++){      
    for(int p=0;p<j;p++) denom[i]=denom[i]+exp(Xbeta[i*j+p]);
    xby[i] = Xbeta[i*j+y[i]-1];
  }
  
  return(sum(xby - log(denom)));
}
