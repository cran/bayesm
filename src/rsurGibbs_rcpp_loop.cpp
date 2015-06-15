#include "bayesm.h"
 
// [[Rcpp::export]]
List rsurGibbs_rcpp_loop(List const& regdata, vec const& indreg, vec const& cumnk, vec const& nk, mat const& XspXs, 
                              mat Sigmainv, mat const& A, vec const& Abetabar, int nu, mat const& V, int nvar, 
                              mat E, mat const& Y, int R, int keep, int nprint){

// Keunwoo Kim 09/19/2014

// Purpose: implement Gibbs Sampler for SUR

// Arguments:
//   Data -- regdata
//           regdata is a list of lists of data for each regression
//           regdata[[i]] contains data for regression equation i
//           regdata[[i]]$y is y, regdata[[i]]$X is X
//           note: each regression can have differing numbers of X vars
//                 but you must have same no of obs in each equation. 
//   Prior -- list of prior hyperparameters
//     betabar,A      prior mean, prior precision
//     nu, V          prior on Sigma
//   Mcmc -- list of MCMC parms
//     R number of draws
//     keep -- thinning parameter
//     nprint - print estimated time remaining on every nprint'th draw

// Output: list of betadraw,Sigmadraw
 
// Model:
//   y_i = X_ibeta + e_i  
//          y is nobs x 1
//          X is nobs x k_i
//          beta is k_i x 1 vector of coefficients
//          i=1,nreg total regressions

//         (e_1,k,...,e_nreg,k) ~ N(0,Sigma) k=1,...,nobs

//   we can also write as stacked regression
//   y = Xbeta+e
//       y is nobs*nreg x 1,X is nobs*nreg x (sum(k_i))
//   routine draws beta -- the stacked vector of all coefficients

// Prior:
//          beta ~ N(betabar,A^-1)
//          Sigma ~ IW(nu,V)

  int reg, mkeep, i, j;
  vec beta, btilde, yti;
  mat IR, ucholinv, EEVinv, Sigma, Xtipyti, Ydti;
  List regdatai, rwout;
  
  int nreg = regdata.size();  
  
  // convert List to std::vector of struct
  std::vector<moments> regdata_vector;
  moments regdatai_struct;
  
  // store vector with struct
  for (reg=0; reg<nreg; reg++){
    regdatai = regdata[reg];
    
    regdatai_struct.y = as<vec>(regdatai["y"]);
    regdatai_struct.X = as<mat>(regdatai["X"]);    
    regdata_vector.push_back(regdatai_struct);    
  }
  
  int nobs = (regdatai_struct.y).size();
  
  mat XtipXti = zeros<mat>(sum(nk), sum(nk));
  mat Sigmadraw(R/keep, nreg*nreg);
  mat betadraw(R/keep, nvar);

  if (nprint>0) startMcmcTimer();

  for (int rep=0; rep<R; rep++){
    
    
    //first draw beta | Sigma
    
    // compute Xtilde'Xtilde
    for (i=0; i<nreg; i++){
      for (j=0; j<nreg; j++){
        XtipXti(span(cumnk[i]-nk[i],cumnk[i]-1), span(cumnk[j]-nk[j],cumnk[j]-1)) =
                  Sigmainv(i,j) * XspXs(span(cumnk[i]-nk[i],cumnk[i]-1), span(cumnk[j]-nk[j],cumnk[j]-1));              
      }      
    }    
    
    // now compute Xtilde'ytilde
    Ydti = Y*Sigmainv;
    Xtipyti = trans(regdata_vector[0].X)*Ydti(span::all,0);
    for (reg=1; reg<nreg; reg++){
      Xtipyti = join_cols(Xtipyti, trans(regdata_vector[reg].X)*Ydti(span::all,reg)); //join_cols is analogous to rbind()
    }     
    
    IR = solve(trimatu(chol(XtipXti + A)), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
    btilde = (IR*trans(IR)) * (Xtipyti + Abetabar);
    beta = btilde + IR*vec(rnorm(nvar));
    
    
    //now draw Sigma | beta
    for (reg=0; reg<nreg; reg++){
      E(span::all,reg) = regdata_vector[reg].y - 
                          regdata_vector[reg].X * beta(span(indreg[reg]-1,indreg[reg+1]-2));
    }
    
    // compute the inverse of E'E+V
    ucholinv = solve(trimatu(chol(trans(E)*E+V)), eye(nreg,nreg));
    EEVinv = ucholinv*trans(ucholinv);
    
    rwout = rwishart(nu+nobs, EEVinv);
    Sigma = as<mat>(rwout["IW"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    Sigmainv = as<mat>(rwout["W"]);
    
    //print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw(mkeep-1, span::all) = trans(beta);
      Sigmadraw(mkeep-1, span::all) = trans(vectorise(Sigma));      
    }
  }
  
  if (nprint>0) endMcmcTimer();
  
  return List::create(
      Named("betadraw") = betadraw,
      Named("Sigmadraw") = Sigmadraw);
}
