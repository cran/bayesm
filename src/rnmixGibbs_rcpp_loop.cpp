#include "bayesm.h"
 
//[[Rcpp::export]]
List rnmixGibbs_rcpp_loop(mat const& y, mat const& Mubar, 
                     mat const& A, int nu, 
                     mat const& V, vec const& a, 
                     vec p, vec z,
                     int const& R, int const& keep, int const& nprint) {

// Wayne Taylor 9/10/2014

  int mkeep = 0;    
  
  mat pdraw(R/keep,p.size());
  mat zdraw(R/keep,z.size());
  List compdraw(R/keep);
  
  if(nprint>0) startMcmcTimer();
  
  // start main iteration loop
  for(int rep = 0; rep<R; rep++) {
    
    List out = rmixGibbs(y, Mubar, A, nu, V, a, p, z);
    
    List compsd = out["comps"];
    p = as<vec>(out["p"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    z = as<vec>(out["z"]);
          
    // print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
            
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      pdraw(mkeep-1,span::all) = trans(p);
      zdraw(mkeep-1,span::all) = trans(z);
      compdraw[mkeep-1] = compsd;
    }
  }
  
  if(nprint>0) endMcmcTimer();
      
  return List::create(
    Named("probdraw") = pdraw, 
    Named("zdraw")    = zdraw, 
    Named("compdraw") = compdraw);
}
