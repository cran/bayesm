#include "bayesm.h"
 
// [[Rcpp::export]]
List rivGibbs_rcpp_loop(vec const& y, vec const& x, mat const& z, mat const& w, vec const& mbg, mat const& Abg, 
                  vec const& md, mat const& Ad, mat const& V, int nu, int R, int keep, int nprint){

// Keunwoo Kim 09/09/2014

// Purpose: draw from posterior for linear I.V. model

// Arguments:
//   Data -- list of z,w,x,y
//        y is vector of obs on lhs var in structural equation
//        x is "endogenous" var in structural eqn
//        w is matrix of obs on "exogenous" vars in the structural eqn
//        z is matrix of obs on instruments
//   Prior -- list of md,Ad,mbg,Abg,nu,V
//        md is prior mean of delta
//        Ad is prior prec
//        mbg is prior mean vector for beta,gamma
//        Abg is prior prec of same
//        nu,V parms for IW on Sigma

//   Mcmc -- list of R,keep 
//        R is number of draws
//        keep is thinning parameter
//        nprint - print estimated time remaining on every nprint'th draw

// Output: list of draws of delta,beta,gamma and Sigma
 
// Model:
//    x=z'delta + e1
//    y=beta*x + w'gamma + e2
//        e1,e2 ~ N(0,Sigma)
//
// Prior:
//   delta ~ N(md,Ad^-1)
//   vec(beta,gamma) ~ N(mbg,Abg^-1)
//   Sigma ~ IW(nu,V)
// 

  vec e1, ee2, bg, u, gamma;
  mat xt, Res, S, B, L, Li, z2, zt1, zt2, ucholinv, VSinv, yt;
  double sig,beta;
  List out;
  int i, mkeep;

  int n = y.size();
  int dimd = z.n_cols;
  int dimg = w.n_cols;

  mat deltadraw(R/keep, dimd);
  vec betadraw(R/keep);
  mat gammadraw(R/keep, dimg);
  mat Sigmadraw(R/keep, 4);  
  mat C = eye(2,2); //eye creates a diagonal matrix

  // set initial values
  mat Sigma = eye(2,2);
  vec delta = 0.1 * ones<vec>(dimd);

  if (nprint>0) startMcmcTimer();  
  
  mat xtd(2*n, dimd);  
  vec zvec = vectorise(trans(z));
  
  // start main iteration loop
  for (int rep=0; rep<R; rep++){   
    
    // draw beta,gamma
    e1 = x - z*delta;
    ee2 = (Sigma(0,1)/Sigma(0,0)) * e1;
    sig = sqrt(Sigma(1,1)-((Sigma(0,1)*Sigma(0,1))/Sigma(0,0)));
    yt = (y-ee2)/sig;
    xt = join_rows(x,w)/sig; //similar to cbind(x,w)
    bg = breg(yt,xt,mbg,Abg);
    beta = bg[0];
    gamma = bg(span(1,bg.size()-1));
    
    // draw delta
    C(1,0) = beta;
    B = C*Sigma*trans(C);
    L = trans(chol(B));
    Li = solve(trimatl(L),eye(2,2)); //trimatl interprets the matrix as lower triangular and makes solve more efficient
    u = y - w*gamma;
    yt = vectorise(Li * trans(join_rows(x,u)));
    z2 = trans(join_rows(zvec, beta*zvec));
    z2 = Li*z2;
    zt1 = z2(0,span::all);
    zt2 = z2(1,span::all);
    zt1.reshape(dimd,n);    
    zt1 = trans(zt1);
    zt2.reshape(dimd,n);    
    zt2 = trans(zt2);
    for (i=0; i<n; i++){
      xtd(2*i,span::all) = zt1(i,span::all);
      xtd(2*i+1,span::all) = zt2(i,span::all);
    }
    delta = breg(yt,xtd,md,Ad);
    
    // draw Sigma
    Res = join_rows(x-z*delta, y-beta*x-w*gamma); //analogous to cbind() 
    S = trans(Res)*Res;
    
    // compute the inverse of V+S
    ucholinv = solve(trimatu(chol(V+S)), eye(2,2));
    VSinv = ucholinv*trans(ucholinv);
    
    out = rwishart(nu+n, VSinv);
    Sigma = as<mat>(out["IW"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    
    // print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      deltadraw(mkeep-1, span::all) = trans(delta);
      betadraw[mkeep-1] = beta;
      gammadraw(mkeep-1, span::all) = trans(gamma);
      Sigmadraw(mkeep-1, span::all) = trans(vectorise(Sigma));
    }    
  }
  
  if (nprint>0) endMcmcTimer();
  
  return List::create(
      Named("deltadraw") = deltadraw,
      Named("betadraw") = NumericVector(betadraw.begin(),betadraw.end()),
      Named("gammadraw") = gammadraw,
      Named("Sigmadraw") = Sigmadraw);   
}




