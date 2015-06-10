#include "bayesm.h"

//[[Rcpp::export]]
vec breg(vec const& y, mat const& X, vec const& betabar, mat const& A) {

// Keunwoo Kim 06/20/2014

// Purpose: draw from posterior for linear regression, sigmasq=1.0

// Output: draw from posterior
 
// Model: y = Xbeta + e  e ~ N(0,I)

// Prior:  beta ~ N(betabar,A^-1)

  int k = betabar.size();
  mat RA = chol(A);
  mat W = join_cols(X, RA); //same as rbind(X,RA)
  vec z = join_cols(y, RA*betabar);
  mat IR = solve(trimatu(chol(trans(W)*W)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  return ((IR*trans(IR))*(trans(W)*z) + IR*vec(rnorm(k)));
} 
