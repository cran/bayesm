#include "bayesm.h"
 
// [[Rcpp::export]]
List rnegbinRw_rcpp_loop(vec const& y, mat const& X, vec const& betabar, mat const& rootA, double a, double b, 
                          vec beta, double alpha, bool fixalpha,
                          mat const& betaroot, double const& alphacroot, int R, int keep, int nprint){

// Keunwoo Kim 11/02/2014

// Arguments:
//       Data
//           X is nobs X nvar matrix
//           y is nobs vector

//       Prior - list containing the prior parameters
//           betabar, rootA - mean of beta prior, chol-root of inverse of variance covariance of beta prior
//           a, b - parameters of alpha prior

//       Mcmc - list containing
//           R is number of draws
//           keep is thinning parameter (def = 1)
//           nprint - print estimated time remaining on every nprint'th draw (def = 100)
//           betaroot - step size for beta RW
//           alphacroot - step size for alpha RW
//           beta - initial guesses for beta
//           alpha - initial guess for alpha
//           fixalpha - if TRUE, fix alpha and draw only beta
//
// Output: 
// 
// Model:
//       (y|lambda,alpha) ~ Negative Binomial(Mean = lambda, Overdispersion par = alpha)
//       ln(lambda) =  X * beta
//
// Prior:
//       beta ~ N(betabar, A^-1)
//       alpha ~ Gamma(a,b) where mean = a/b and variance = a/(b^2)
//
  vec betac;
  double ldiff, acc, unif, logalphac, oldlpostalpha, oldlpostbeta, clpostbeta, clpostalpha;
  int mkeep, rep;
  
  int nvar = X.n_cols;  
  int nacceptbeta = 0;
  int nacceptalpha = 0;  

  vec alphadraw(R/keep);
  mat betadraw(R/keep, nvar);
  
  if (nprint>0) startMcmcTimer();
  
  //start main iteration loop
  for (rep=0; rep<R; rep++){
    
    // Draw beta
    betac = beta + betaroot*vec(rnorm(nvar));
    oldlpostbeta = lpostbeta(alpha, beta, X, y, betabar, rootA);
    clpostbeta = lpostbeta(alpha, betac, X, y, betabar, rootA);
    ldiff = clpostbeta - oldlpostbeta;
    acc = exp(ldiff);
    if (acc > 1) acc = 1;    
    if(acc < 1) {unif=runif(1)[0];} else {unif=0;} //runif returns a NumericVector, so using [0] allows for conversion to double by extracting the first element
    if (unif <= acc){
      beta = betac;
      nacceptbeta = nacceptbeta + 1;
    } 
    
    // Draw alpha
    if (!fixalpha){
      logalphac = log(alpha) + alphacroot*rnorm(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
      oldlpostalpha = lpostalpha(alpha, beta, X, y, a, b);
      clpostalpha = lpostalpha(exp(logalphac), beta, X, y, a, b);
      ldiff = clpostalpha - oldlpostalpha;
      acc = exp(ldiff);
      if (acc > 1) acc = 1;    
      if(acc < 1) {unif=runif(1)[0];} else {unif=0;} //runif returns a NumericVector, so using [0] allows for conversion to double by extracting the first element
      if (unif <= acc){
        alpha = exp(logalphac);
        nacceptalpha = nacceptalpha + 1;
      }
    }

    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);    
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw(mkeep-1, span::all) = trans(beta);
      alphadraw[mkeep-1] = alpha;           
    } 
  }
    
  if (nprint>0) endMcmcTimer();
  return List::create(
      Named("betadraw") = betadraw,
      Named("alphadraw") = alphadraw,      
      Named("nacceptbeta") = nacceptbeta,
      Named("nacceptalpha") = nacceptalpha);
}
