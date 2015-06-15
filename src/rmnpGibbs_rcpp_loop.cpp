#include "bayesm.h"
 
//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------
vec drawwi(vec const& w, vec const& mu, mat const& sigmai, int p, int y){

// Wayne Taylor 9/8/2014

//function to draw w_i by Gibbing thru p vector

  int above;
	double bound;
  vec outwi = w;
  vec maxInd(2);

	for(int i = 0; i<p; i++){	
		bound = 0.0;
		for(int j = 0; j<p; j++) if(j!=i) {
        maxInd[0] = bound;
        maxInd[1] = outwi[j];
        bound = max(maxInd);}
    
    if (y==(i+1))
			above = 0;
		else 
			above = 1;
    
		vec CMout = condmom(outwi,mu,sigmai,p,i+1);
    outwi[i] = rtrun1(CMout[0],CMout[1],bound,above);
  }

  return (outwi);
}

vec draww(vec const& w, vec const& mu, mat const& sigmai, ivec const& y){

// Wayne Taylor 9/8/2014 

//function to gibbs down entire w vector for all n obs
  
  int n = y.n_rows;
  int p = sigmai.n_cols;
  int ind; 
  vec outw = zeros<vec>(w.n_rows);
  
	for(int i = 0; i<n; i++){
    ind = p*i;
		outw.subvec(ind,ind+p-1) = drawwi(w.subvec(ind,ind+p-1),mu.subvec(ind,ind+p-1),sigmai,p,y[i]);
	}

  return (outw);
}

//MAIN FUNCTION---------------------------------------------------------------------------------------
//[[Rcpp::export]]
List rmnpGibbs_rcpp_loop(int R, int keep, int nprint, int pm1, 
                         ivec const& y, mat const& X, vec const& beta0, mat const& sigma0, 
                         mat const& V, int nu, vec const& betabar, mat const& A) {

// Wayne Taylor 9/24/2014

  int n = y.n_rows;
  int k = X.n_cols;
  int Xrows = X.n_rows;
  
  //allocate space for draws
  mat sigmadraw = zeros<mat>(R/keep, pm1*pm1);
  mat betadraw = zeros<mat>(R/keep,k);
  vec wnew = zeros<vec>(Xrows);
  
  //set initial values of w,beta, sigma (or root of inv)
  vec wold = wnew;
  vec betaold = beta0;
  
  mat C = chol(solve(trimatu(sigma0),eye(sigma0.n_cols,sigma0.n_cols))); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  //C is upper triangular root of sigma^-1 (G) = C'C
  
  mat sigmai, zmat, epsilon, S, IW, ucholinv, VSinv;
  vec betanew;
  List W;
  
  // start main iteration loop
  int mkeep = 0;
  
  if(nprint>0) startMcmcTimer();
  
    for(int rep = 0; rep<R; rep++) {
      
      //draw w given beta(rep-1),sigma(rep-1)
      sigmai = trans(C)*C;
      //    draw latent vector
    
      //    w is n x (p-1) vector
      //       X ix n(p-1) x k  matrix
      //       y is multinomial 1,..., p
      //       beta is k x 1 vector
      //       sigmai is (p-1) x (p-1) 
          
      wnew = draww(wold,X*betaold,sigmai,y);
      
      //draw beta given w(rep) and sigma(rep-1)
      //  note:  if Sigma^-1 (G) = C'C then Var(Ce)=CSigmaC' = I
      //  first, transform w_i = X_ibeta + e_i by premultiply by C
      
      zmat = join_rows(wnew,X);
      zmat.reshape(pm1,n*(k+1));
      zmat = C*zmat;
      zmat.reshape(Xrows,k+1);
      
      betanew = breg(zmat(span::all,0),zmat(span::all,span(1,k)),betabar,A);
      
      //draw sigmai given w and beta
      epsilon = wnew-X*betanew;
      epsilon.reshape(pm1,n);  
      S = epsilon*trans(epsilon);
      
      //same as chol2inv(chol(V+S))
      ucholinv = solve(trimatu(chol(V+S)), eye(S.n_cols,S.n_cols));
      VSinv = ucholinv*trans(ucholinv);
      
      W = rwishart(nu+n,VSinv);
      C = as<mat>(W["C"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      
      //print time to completion
      if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
      
      //save every keepth draw
        if((rep+1)%keep==0){
          mkeep = (rep+1)/keep;
          betadraw(mkeep-1,span::all) = trans(betanew);
          IW  = as<mat>(W["IW"]);
          sigmadraw(mkeep-1,span::all) = trans(vectorise(IW));
         }
        
      wold = wnew;
      betaold = betanew;
    }
  
  if(nprint>0) endMcmcTimer();
      
  return List::create(
    Named("betadraw") = betadraw, 
    Named("sigmadraw") = sigmadraw);
}
