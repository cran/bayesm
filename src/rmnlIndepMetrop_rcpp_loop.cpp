#include "bayesm.h"
 
//[[Rcpp::export]]
List rmnlIndepMetrop_rcpp_loop(int R, int keep, int nu,
                                vec const& betastar, mat const& root,vec const& y,mat const& X,
                                vec const& betabar,mat const& rootpi,mat const& rooti,
                                double oldlimp,double oldlpost,int nprint) {

// Wayne Taylor 9/7/2014

  int mkeep = 0;
  int naccept = 0;    
  int ncolX = X.n_cols;
  
  mat betadraw(R/keep, ncolX);
  vec loglike(R/keep);
  vec betac = zeros<vec>(ncolX);
  rowvec beta = zeros<rowvec>(ncolX);
  double cloglike, clpost, climp, ldiff, alpha, unif, oldloglike;
  vec alphaminv;
  
  if(nprint>0) startMcmcTimer();
  
  // start main iteration loop
  for(int rep = 0; rep<R; rep++) {
    
    betac = rmvst(nu,betastar,root);
    cloglike = llmnl(betac,y,X);
    clpost = cloglike+lndMvn(betac,betabar,rootpi);
    climp = lndMvst(betac,nu,betastar,rooti,false);
    ldiff = clpost+oldlimp-oldlpost-climp;
    alphaminv << 1 << exp(ldiff); //intializes variables in the alphaminv vec: c(1,exp(ldiff))
    alpha = min(alphaminv);
  
    if(alpha < 1.0) {
        unif = runif(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
      }else{
        unif = 0.0;
      }
    if (unif <= alpha){ 
      beta = trans(betac);
      oldloglike = cloglike;
      oldlpost = clpost;
      oldlimp = climp;
      naccept++;
    }
          
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw(mkeep-1,span::all) = beta;
      loglike[mkeep-1] = oldloglike;
    }
  }
  
  if(nprint>0) endMcmcTimer();
      
  return List::create(
    Named("betadraw") = betadraw, 
    Named("loglike") = loglike, 
    Named("naccept") = naccept);
}
