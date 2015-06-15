#include "bayesm.h"
 
//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------
//dstartoc is a fuction to transfer dstar to its cut-off value    
vec dstartoc(vec const& dstar){
  int ndstar = dstar.size();
  vec c(ndstar+3);
  c[0] = -100;
  c[1] = 0;
  c(span(2,ndstar+1)) = cumsum(exp(dstar));
  c[ndstar+2] = 100;
  
  return (c);
} 

// compute conditional likelihood of data given cut-offs
double lldstar(vec const& dstar, vec const& y, vec const& mu){
  vec gamma = dstartoc(dstar);
  
  int ny = y.size();
  NumericVector gamma1(ny);
  NumericVector gamma2(ny);
  for (int i=0; i<ny; i++){
    gamma1[i] = gamma(y[i]);
    gamma2[i] = gamma(y[i]-1);
  }
  NumericVector temp = pnorm(gamma1-as<NumericVector>(wrap(mu)))-pnorm(gamma2-as<NumericVector>(wrap(mu))); //pnorm takes Rcpp type NumericVector, NOT arma objects of type vec
  vec arg = as<vec>(temp);
  double epsilon = 1.0/(10^-50);
  for (int j=0; j<ny; j++){
    if (arg[j]<epsilon){
      arg[j] = epsilon;
    }
  }
  return (sum(log(arg)));
}

List dstarRwMetrop(vec const& y, vec const& mu, vec const& olddstar, double s, mat const& inc_root, 
                    vec const& dstarbar, double oldll, mat const& rootdi, int ncut){ 

// function to execute rw metropolis for the dstar
// y is n vector with element = 1,...,j 
// X is n x k matrix of x values 
// RW increments are N(0,s^2*t(inc.root)%*%inc.root)
// prior on dstar is N(dstarbar,Sigma)  Sigma^-1=rootdi*t(rootdi)
//  inc.root, rootdi are upper triangular
//  this means that we are using the UL decomp of Sigma^-1 for prior 
// olddstar is the current
//
  int stay = 0;
  double unif;
  vec dstardraw;

  vec dstarc = olddstar + s*trans(inc_root)*vec(rnorm(ncut));
  double cll = lldstar(dstarc, y, mu);
  double clpost = cll + lndMvn(dstarc, dstarbar, rootdi);
  double ldiff = clpost - oldll - lndMvn(olddstar, dstarbar, rootdi);
  double alpha = exp(ldiff);
  
  if (alpha>1){
    alpha = 1.0;
  } 

  if (alpha<1){
    unif = runif(1)[0]; //runif returns a NumericVector, so using [0] allows for conversion to double by extracting the first element
  }
  else{
    unif = 0;
  }
  
  if (unif<=alpha){
    dstardraw = dstarc; 
    oldll = cll;
  }
  else{
    dstardraw = olddstar;
    stay = 1;
  }
  
  return List::create(
      Named("dstardraw") = dstardraw,
      Named("oldll") = oldll,
      Named("stay") = stay
  );
}   

//MAIN FUNCTION---------------------------------------------------------------------------------------
// [[Rcpp::export]]
List rordprobitGibbs_rcpp_loop(vec const& y, mat const& X, int k, mat const& A, vec const& betabar, mat const& Ad, 
                          double s, mat const& inc_root, vec const& dstarbar, vec const& betahat, 
                          int R, int keep, int nprint){

// Keunwoo Kim 09/09/2014

// Purpose: draw from posterior for ordered probit using Gibbs Sampler and metropolis RW

// Arguments:
//  Data
//    X is nobs x nvar, y is nobs vector of 1,2,.,k (ordinal variable)
//  Prior
//    A is nvar x nvar prior preci matrix
//    betabar is nvar x 1 prior mean
//    Ad is ndstar x ndstar prior preci matrix of dstar (ncut is number of cut-offs being estimated)
//    dstarbar is ndstar x 1 prior mean of dstar
//  Mcmc
//    R is number of draws
//    keep is thinning parameter
//    nprint - prints the estimated time remaining for every nprint'th draw
//    s is scale parameter of random work Metropolis

// Output: list of betadraws and cutdraws
 
// Model: 
//    z=Xbeta + e  < 0  e ~N(0,1)
//    y=1,..,k, if z~c(c[k], c[k+1])

//    cutoffs = c[1],..,c[k+1]
//    dstar = dstar[1],dstar[k-2]
//    set c[1]=-100, c[2]=0, ...,c[k+1]=100

//    c[3]=exp(dstar[1]),c[4]=c[3]+exp(dstar[2]),...,
//    c[k]=c[k-1]+exp(datsr[k-2])
    
// Note: 1. length of dstar = length of cutoffs - 3
//       2. Be careful in assessing prior parameter, Ad.  .1 is too small for many applications.

// Prior: 
//  beta ~ N(betabar,A^-1)
//  dstar ~ N(dstarbar, Ad^-1)

  int stay, i, mkeep;
  vec z;
  List metropout;
 
  int nvar = X.n_cols;
  int ncuts = k+1;
  int ncut = ncuts-3;
  int ndstar = k-2;
  int ny = y.size();

  mat betadraw(R/keep, nvar);
  mat cutdraw(R/keep, ncuts);
  mat dstardraw(R/keep, ndstar);
  vec staydraw(R/keep);
  vec cutoff1(ny);
  vec cutoff2(ny);
  vec sigma(X.n_rows); sigma.ones();
  
  // compute the inverse of trans(X)*X+A
  mat ucholinv = solve(trimatu(chol(trans(X)*X+A)), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  mat XXAinv = ucholinv*trans(ucholinv);

  mat root = chol(XXAinv);
  vec Abetabar = trans(A)*betabar;
  
  // compute the inverse of Ad
  ucholinv = solve(trimatu(chol(Ad)), eye(ndstar,ndstar));
  mat Adinv = ucholinv*trans(ucholinv);
  
  mat rootdi = chol(Adinv);
  
  // set initial values for MCMC  
  vec olddstar(ndstar); 
  olddstar.zeros();
  vec beta = betahat;    
  vec cutoffs = dstartoc(olddstar);  
  double oldll = lldstar(olddstar, y, X*betahat);
  
  if (nprint>0) startMcmcTimer();
  
  //start main iteration loop
  for (int rep=0; rep<R; rep++){
    
    //draw z given beta(i-1), sigma, y, cut-offs
    for (i=0; i<ny; i++){
      cutoff1[i] = cutoffs[y[i]-1];
      cutoff2[i] = cutoffs[y[i]];
    }
    z = rtrunVec(X*beta, sigma, cutoff1, cutoff2);

    //draw beta given z and rest
    beta = breg1(root,X,z,Abetabar);
   
    //draw gamma given z
    metropout = dstarRwMetrop(y,X*beta,olddstar,s,inc_root,dstarbar,oldll,rootdi, ncut);   
    olddstar = as<vec>(metropout["dstardraw"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    oldll =  as<double>(metropout["oldll"]);
    cutoffs = dstartoc(olddstar);
    stay = as<int>(metropout["stay"]);  

    //print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      cutdraw(mkeep-1,span::all) = trans(cutoffs);
      dstardraw(mkeep-1,span::all) = trans(olddstar);
      betadraw(mkeep-1,span::all) = trans(beta);
      staydraw[mkeep-1] = stay;
    }                
  }
  double accept = 1-sum(staydraw)/(R/keep);
  if (nprint>0) endMcmcTimer();

  return List::create(
      Named("cutdraw") = cutdraw,
      Named("dstardraw") = dstardraw,
      Named("betadraw") = betadraw,
      Named("accept") = accept
  );
}
