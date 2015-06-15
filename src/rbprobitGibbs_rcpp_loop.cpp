#include "bayesm.h"
 
// [[Rcpp::export]]
List rbprobitGibbs_rcpp_loop(vec const& y, mat const& X, vec const& Abetabar, mat const& root, 
                        vec beta, vec const& sigma, vec const& a, vec const& b, int R, int keep, int nprint){

// Keunwoo Kim 09/09/2014

// Purpose: draw from posterior for binary probit using Gibbs Sampler

// Arguments:
//  X is nobs x nvar, y is nobs vector of 0,1
//  A is nvar x nvar prior preci matrix
//  betabar is nvar x 1 prior mean
//  R is number of draws
//  keep is thinning parameter
//  nprint - prints the estimated time remaining for every nprint'th draw

// Output: list of betadraws
 
// Model: y = 1 if  w=Xbeta+e>0  e~N(0,1)

// Prior: beta ~ N(betabar,A^-1)
 
  int mkeep;
  vec mu;
  vec z;

  int nvar = X.n_cols;
  
  mat betadraw(R/keep, nvar);
  
  if (nprint>0) startMcmcTimer();
  
  //start main iteration loop
  for (int rep=0; rep<R; rep++){
    
    // draw z given beta(i-1)
    mu = X*beta;
    z = rtrunVec(mu, sigma, a, b);
    beta = breg1(root, X, z, Abetabar);

    // print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw(mkeep-1, span::all) = trans(beta);      
    }
  }
  
  if (nprint>0) endMcmcTimer();
  
  return List::create(Named("betadraw") = betadraw);
}
