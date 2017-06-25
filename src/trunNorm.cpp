// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <float.h>
#include <Rmath.h>
#include <math.h>
#include "bayesm.h"
using namespace Rcpp;
using namespace arma;


double dexpr(double const& a) {
  // rossi 10/2016
  // routine to draw use rejection sampling to draw from a truncated normal
  // truncated to the tail.  Use exponential distribution as the acceptance function
  //  a is the rejection point (usually, (a-mu)/sigma)) 
  double x,e,e1; 
  int success;
  success=0;
  while(success == 0){
    e = -log(runif(1)[0]);
    e1 = -log(runif(1)[0]);
    if(pow(e,2) <= 2.0*e1*pow(a,2)){
      x = a + e/a;
      success = 1;
    }
  }
  return(x);
}

double invCdfNorm(double const& a) {
  // rossi 10/2016
  // draw from truncated normal truncated from below by a using inverse CDF method
  double Phia,unif,arg,z;
  Phia = R::pnorm(a,0.0,1.0,1,0);
  unif = runif(1)[0];
  arg = unif*(1.0-Phia)+Phia;
  z = R::qnorm(arg,0.0,1.0,1,0);
  return(z);
}


double dnr(double const& a) {
  // rossi 10/2016
  // draw from standard normal truncated below by a using rejection sampling 
  double candz,z; 
  int success;
  success=0;
  while(success == 0){
    candz=rnorm(1)[0];
    if(candz >= a){
      z=candz; 
      success =1;
    }
    }
  return(z);
}

double trunNormBelow(double const& a){
  // rossi 10/2016
  // draw from standard normal truncated below by a
  // we divide into three regions (regions given for trun below)
  // let a = (trunpt-mu)/sig
  // 1: a > 4, use exponential rejection sampling
  // 2: -4 < a <=4, use inverse CDF
  // 3: a <= -4, normal rejection sampling
  double z;
  if (a > 4){
    // do tail sampling using exponential rejection
    z=dexpr(a);
  }
  else {
    if (a <= -4){
      // normal rejection sampling
      z=dnr(a);
    }
    // -4 < a <=4
    z=invCdfNorm(a); 
  }
 return(z); 
}

double trunNorm(double mu,double sig, double trunpt, int above){
  // rossi 10/2016
  // function to draw from N(mu,sig)  truncated by trunpt
  // if above=0, truncate from below at trunpt
  // if above=1, truncate from above at trunpt
  double a,z,draw;
  if(above == 0){
    a=(trunpt-mu)/sig;
    z= trunNormBelow(a);
    draw = sig*z + mu;
  }
  else {
    a=(mu-trunpt)/sig;
    z= trunNormBelow(a);
    draw = -sig*z + mu;
  }
  return(draw);
}


// vec trunNorm_multi(double mu, double sig, double trunpt, int above, int nd){
//    // Dan Yavorsky 10/16
//    // Just loops over trunNorm, multiple draws from single truncated distribution
//    vec rtn_vec(nd);
//    for (int i = 0; i<nd; i++) {
//        rtn_vec[i] = trunNorm(mu, sig, trunpt, above);
//    }
//    return(rtn_vec);
// }

vec trunNorm_vec(vec const& mu, vec const& sig, vec const& trunpt, vec const& above){
  // Dan Yavorsky 10/16
  // Loops over trunNorm, everything is vectorized
   int nd = mu.size();
   vec rtn_vec(nd);
   for (int i = 0; i<nd; i++) {
       rtn_vec[i] = trunNorm(mu[i], sig[i], trunpt[i], above[i]);
   }
   return(rtn_vec);
}






