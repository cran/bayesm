#include "bayesm.h"
 
//[[Rcpp::export]]
List rhierLinearMixture_rcpp_loop(List const& regdata, mat const& Z,
                                  vec const& deltabar, mat const& Ad, mat const& mubar, mat const& Amu,
                                  int const& nu, mat const& V, int nu_e, vec const& ssq,
                                  int R, int keep, int nprint, bool drawdelta,
                                  mat olddelta,  vec const& a, vec oldprob, vec ind, vec tau){

// Wayne Taylor 10/02/2014

  int nreg = regdata.size();
  int nvar = V.n_cols;
  int nz = Z.n_cols;
  
  mat rootpi, betabar, Abeta, Abetabar;
  int mkeep;
  unireg runiregout_struct;
  List regdatai, nmix;
  
  // convert List to std::vector of type "moments"
  std::vector<moments> regdata_vector;
  moments regdatai_struct;
  
  // store vector with struct
  for (int reg = 0; reg<nreg; reg++){
    regdatai = regdata[reg];
    
    regdatai_struct.y = as<vec>(regdatai["y"]);
    regdatai_struct.X = as<mat>(regdatai["X"]);
    regdatai_struct.XpX = as<mat>(regdatai["XpX"]);
    regdatai_struct.Xpy = as<vec>(regdatai["Xpy"]);
    regdata_vector.push_back(regdatai_struct);    
  }
  
  // allocate space for draws
  mat oldbetas = zeros<mat>(nreg,nvar);
  mat taudraw(R/keep, nreg);
  cube betadraw(nreg, nvar, R/keep);
  mat probdraw(R/keep, oldprob.size());
  mat Deltadraw(1,1); if(drawdelta) Deltadraw.zeros(R/keep, nz*nvar);//enlarge Deltadraw only if the space is required
  List compdraw(R/keep);
  
  if (nprint>0) startMcmcTimer();

  for (int rep = 0; rep<R; rep++){
   
   //first draw comps,ind,p | {beta_i}, delta
   // ind,p need initialization comps is drawn first in sub-Gibbs
   List mgout;
   if(drawdelta) {
      olddelta.reshape(nvar,nz);
      mgout = rmixGibbs(oldbetas-Z*trans(olddelta),mubar,Amu,nu,V,a,oldprob,ind);
    } else {
      mgout = rmixGibbs(oldbetas,mubar,Amu,nu,V,a,oldprob,ind);
    }
   
   List oldcomp = mgout["comps"];
   oldprob = as<vec>(mgout["p"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
   ind = as<vec>(mgout["z"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
   
  //now draw delta | {beta_i}, ind, comps
   if(drawdelta) olddelta = drawDelta(Z,oldbetas,ind,oldcomp,deltabar,Ad);
   
  //loop over all regression equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
      for(int reg = 0; reg<nreg; reg++){
        List oldcompreg = oldcomp[ind[reg]-1];
        rootpi = as<mat>(oldcompreg[1]);
        
        //note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
        if(drawdelta){
          olddelta.reshape(nvar,nz);
          betabar = as<vec>(oldcompreg[0])+olddelta*vectorise(Z(reg,span::all));
        } else {
          betabar = as<vec>(oldcompreg[0]);
        }
      
        Abeta = trans(rootpi)*rootpi;
        Abetabar = Abeta*betabar;

        runiregout_struct = runiregG(regdata_vector[reg].y, regdata_vector[reg].X,
                                regdata_vector[reg].XpX, regdata_vector[reg].Xpy, 
                                tau[reg], Abeta, Abetabar, nu_e, ssq[reg]);
      
        oldbetas(reg,span::all) = trans(runiregout_struct.beta);
        tau[reg] = runiregout_struct.sigmasq;
      }
      
  //print time to completion and draw # every nprint'th draw
  if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      taudraw(mkeep-1, span::all) = trans(tau);
      betadraw.slice(mkeep-1) = oldbetas;
      probdraw(mkeep-1, span::all) = trans(oldprob);
      if(drawdelta) Deltadraw(mkeep-1, span::all) = trans(vectorise(olddelta));
      compdraw[mkeep-1] = oldcomp;
    }
  }
  
  if (nprint>0) endMcmcTimer();
  
  nmix = List::create(Named("probdraw") = probdraw,
  				  Named("zdraw") = R_NilValue, //sets the value to NULL in R
					  Named("compdraw") = compdraw);
	
  if(drawdelta){
    return(List::create(
      Named("taudraw") = taudraw,
      Named("Deltadraw") = Deltadraw,
      Named("betadraw") = betadraw,
      Named("nmix") = nmix));
	} else {
    return(List::create(
      Named("taudraw") = taudraw,
      Named("betadraw") = betadraw,
      Named("nmix") = nmix));
  }
}
