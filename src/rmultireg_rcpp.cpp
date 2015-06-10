#include "bayesm.h"
 
// [[Rcpp::export]]
List rmultireg(mat const& Y, mat const& X, mat const& Bbar, mat const& A, int nu, mat const& V) {

// Keunwoo Kim 09/09/2014

// Purpose: draw from posterior for Multivariate Regression Model with natural conjugate prior

// Arguments:
//  Y is n x m matrix
//  X is n x k
//  Bbar is the prior mean of regression coefficients  (k x m)
//  A is prior precision matrix
//  nu, V are parameters for prior on Sigma

// Output: list of B, Sigma draws of matrix of coefficients and Sigma matrix
 
// Model: 
//  Y=XB+U  cov(u_i) = Sigma
//  B is k x m matrix of coefficients

// Prior:  
//  beta|Sigma  ~ N(betabar,Sigma (x) A^-1)
//  betabar=vec(Bbar)
//  beta = vec(B) 
//  Sigma ~ IW(nu,V) or Sigma^-1 ~ W(nu, V^-1)

  int n = Y.n_rows;
  int m = Y.n_cols;
  int k = X.n_cols;
  
  //first draw Sigma
  mat RA = chol(A);
  mat W = join_cols(X, RA); //analogous to rbind() in R
  mat Z = join_cols(Y, RA*Bbar);
  // note:  Y,X,A,Bbar must be matrices!
  mat IR = solve(trimatu(chol(trans(W)*W)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  // W'W = R'R  &  (W'W)^-1 = IRIR'  -- this is the UL decomp!
  mat Btilde = (IR*trans(IR)) * (trans(W)*Z);
  // IRIR'(W'Z) = (X'X+A)^-1(X'Y + ABbar)
  mat E = Z-W*Btilde;
  mat S = trans(E)*E;
  // E'E
  
  // compute the inverse of V+S
  mat ucholinv = solve(trimatu(chol(V+S)), eye(m,m));
  mat VSinv = ucholinv*trans(ucholinv);
  
  List rwout = rwishart(nu+n, VSinv);
  
  // now draw B given Sigma
  //   note beta ~ N(vec(Btilde),Sigma (x) Covxxa)
  //       Cov=(X'X + A)^-1  = IR t(IR)  
  //       Sigma=CICI'    
  //       therefore, cov(beta)= Omega = CICI' (x) IR IR' = (CI (x) IR) (CI (x) IR)'
  //  so to draw beta we do beta= vec(Btilde) +(CI (x) IR)vec(Z_mk)  
  //  		Z_mk is m x k matrix of N(0,1)
  //	since vec(ABC) = (C' (x) A)vec(B), we have 
  //		B = Btilde + IR Z_mk CI'

  mat CI = rwout["CI"]; //there is no need to use as<mat>(rwout["CI"]) since CI is being initiated as a mat in the same line
  mat draw = mat(rnorm(k*m));
  draw.reshape(k,m);
  mat B = Btilde + IR*draw*trans(CI);
    
  return List::create(
      Named("B") = B, 
      Named("Sigma") = rwout["IW"]);
}
