#include "bayesm.h"

// [[Rcpp::export]]
List rhierLinearModel_rcpp_loop(List const& regdata, mat const& Z, mat const& Deltabar, mat const& A, int nu, 
                          mat const& V, int nu_e, vec const& ssq, vec tau, mat Delta, mat Vbeta, int R, int keep, int nprint){

// Keunwoo Kim 09/16/2014

// Purpose: run hiearchical regression model

// Arguments:
//   Data list of regdata,Z 
//     regdata is a list of lists each list with members y, X
//        e.g. regdata[[i]]=list(y=y,X=X)
//     X has nvar columns
//     Z is nreg=length(regdata) x nz

//   Prior list of prior hyperparameters
//     Deltabar,A, nu.e,ssq,nu,V
//          note: ssq is a nreg x 1 vector!

//   Mcmc
//     list of Mcmc parameters
//     R is number of draws
//     keep is thining parameter -- keep every keepth draw
//     nprint - print estimated time remaining on every nprint'th draw

// Output: 
//   list of 
//   betadraw -- nreg x nvar x R/keep array of individual regression betas
//   taudraw -- R/keep x nreg  array of error variances for each regression
//   Deltadraw -- R/keep x nz x nvar array of Delta draws
//   Vbetadraw -- R/keep x nvar*nvar array of Vbeta draws

// Model:
// nreg regression equations 
//        y_i = X_ibeta_i + epsilon_i  
//        epsilon_i ~ N(0,tau_i)
//             nvar X vars in each equation

// Prior:
//        tau_i ~ nu.e*ssq_i/chisq(nu.e)  tau_i is the variance of epsilon_i
//        beta_i ~ N(ZDelta[i,],V_beta)
//               Note:  ZDelta is the matrix Z * Delta; [i,] refers to ith row of this product!

//          vec(Delta) | V_beta ~ N(vec(Deltabar),Vbeta (x) A^-1)
//          V_beta ~ IW(nu,V)  or V_beta^-1 ~ W(nu,V^-1)
//              Delta, Deltabar are nz x nvar
//              A is nz x nz
//              Vbeta is nvar x nvar
        
//          NOTE: if you don't have any z vars, set Z=iota (nreg x 1)
 
// Update Note:
//        (Keunwoo Kim 04/07/2015)
//        Changed "rmultireg" to return List object, which is the original function.
//        Efficiency is almost same as when the output is a struct object.
//        Nothing different from "rmultireg1" in the previous R version.

  int reg, mkeep;
  mat Abeta, betabar, ucholinv, Abetabar;
  List regdatai, rmregout;
  unireg regout_struct;
  
  int nreg = regdata.size();
  int nvar = V.n_cols;
  int nz = Z.n_cols;
  
  // convert List to std::vector of struct
  std::vector<moments> regdata_vector;
  moments regdatai_struct;
  
  // store vector with struct
  for (reg=0; reg<nreg; reg++){
    regdatai = regdata[reg];
    
    regdatai_struct.y = as<vec>(regdatai["y"]);
    regdatai_struct.X = as<mat>(regdatai["X"]);
    regdatai_struct.XpX = as<mat>(regdatai["XpX"]);
    regdatai_struct.Xpy = as<vec>(regdatai["Xpy"]);
    regdata_vector.push_back(regdatai_struct);    
  } 
  
  mat betas(nreg, nvar);
  mat Vbetadraw(R/keep, nvar*nvar);
  mat Deltadraw(R/keep, nz*nvar);
  mat taudraw(R/keep, nreg);
  cube betadraw(nreg, nvar, R/keep);

  if (nprint>0) startMcmcTimer();
  
  //start main iteration loop
  for (int rep=0; rep<R; rep++){    

    // compute the inverse of Vbeta
    ucholinv = solve(trimatu(chol(Vbeta)), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
    Abeta = ucholinv*trans(ucholinv);
    
    betabar = Z*Delta;
    Abetabar = Abeta*trans(betabar);
    
    //loop over all regressions
    for (reg=0; reg<nreg; reg++){      
    
      regout_struct = runiregG(regdata_vector[reg].y, regdata_vector[reg].X, 
                                regdata_vector[reg].XpX, regdata_vector[reg].Xpy, 
                                tau[reg], Abeta, Abetabar(span::all,reg), 
                                nu_e, ssq[reg]);
      betas(reg,span::all) = trans(regout_struct.beta);
      tau[reg] = regout_struct.sigmasq;
    }
    
    //draw Vbeta, Delta | {beta_i}
    rmregout = rmultireg(betas,Z,Deltabar,A,nu,V);
    Vbeta = as<mat>(rmregout["Sigma"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    Delta = as<mat>(rmregout["B"]);
  
    //print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      Vbetadraw(mkeep-1, span::all) = trans(vectorise(Vbeta));
      Deltadraw(mkeep-1, span::all) = trans(vectorise(Delta));
      taudraw(mkeep-1, span::all) = trans(tau);
      betadraw.slice(mkeep-1) = betas;
    }    
  }
  
  if (nprint>0) endMcmcTimer();
  
  return List::create(
    Named("Vbetadraw") = Vbetadraw,
    Named("Deltadraw") = Deltadraw,
	  Named("betadraw") = betadraw,
	  Named("taudraw") = taudraw);
}
