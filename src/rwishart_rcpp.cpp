#include "bayesm.h"
 
// [[Rcpp::export]]
List rwishart(int const& nu, mat const& V){

// Wayne Taylor 4/7/2015

// Function to draw from Wishart (nu,V) and IW
 
// W ~ W(nu,V)
// E[W]=nuV

// WI=W^-1
// E[WI]=V^-1/(nu-m-1)
  
  // T has sqrt chisqs on diagonal and normals below diagonal
  int m = V.n_rows;
  mat T = zeros(m,m);
  
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]); //rchisq returns a vectorized object, so using [0] allows for the conversion to double
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = rnorm(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
  }}
  
  mat C = trans(T)*chol(V);
  mat CI = solve(trimatu(C),eye(m,m)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  // C is the upper triangular root of Wishart therefore, W=C'C
  // this is the LU decomposition Inv(W) = CICI' Note: this is
  // the UL decomp not LU!
  
  // W is Wishart draw, IW is W^-1
  
  return List::create(
    Named("W") = trans(C) * C,
    Named("IW") = CI * trans(CI),
    Named("C") = C,
    Named("CI") = CI);
}
