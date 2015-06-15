#include "bayesm.h"
 
//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------
double llnegbinpooled(std::vector<moments> regdata_vector, mat Beta, double alpha){
  
// Wayne Taylor 12/01/2014

// "Unlists" the regdata and calculates the negative binomial loglikelihood using individual-level betas
  
  int nreg = regdata_vector.size();
  double ll = 0.0;
  
  for(int reg = 0; reg<nreg; reg++){
  vec lambda = exp(regdata_vector[reg].X*trans(Beta(reg,span::all)));
  ll = ll + llnegbin(regdata_vector[reg].y,lambda,alpha,TRUE);
  }
  
  return(ll);
}

// [[Rcpp::export]]
List rhierNegbinRw_rcpp_loop(List const& regdata, List const& hessdata, mat const& Z, mat Beta, mat Delta,
                             mat const& Deltabar, mat const& Adelta, int nu, mat const& V, double a, double b,
                             int R, int keep, double sbeta, double alphacroot, int nprint, mat rootA,
                             double alpha, bool fixalpha){
                            
// Wayne Taylor 12/01/2014                          

//   Model
//       (y_i|lambda_i,alpha) ~ Negative Binomial(Mean = lambda_i, Overdispersion par = alpha)
//
//       ln(lambda_i) =  X_i * beta_i
//
//       beta_i = Delta'*z_i + nu_i
//               nu_i~N(0,Vbeta)
//       Note: rootA = the Cholesky root of the inverse of Vbeta
//
//   Priors
//       vec(Delta|Vbeta) ~ N(vec(Deltabar), Vbeta (x) (Adelta^-1))
//       Vbeta ~ Inv Wishart(nu, V)
//       alpha ~ Gamma(a,b) where mean = a/b and variance = a/(b^2)
//
//   Arguments
//       Data = list of regdata,Z 
//           regdata is a list of lists each list with members y, X
//              e.g. regdata[[i]]=list(y=y,X=X)
//              X has nvar columns including a first column of ones
//              Z is nreg=length(regdata) x nz with a first column of ones
//
//       Prior - list containing the prior parameters
//           Deltabar, Adelta - mean of Delta prior, inverse of variance covariance of Delta prior
//           nu, V - parameters of Vbeta prior
//           a, b - parameters of alpha prior
//
//       Mcmc - list containing
//           R is number of draws
//           keep is thinning parameter (def = 1)
//           nprint - print estimated time remaining on every nprint'th draw (def = 100)
//           s_beta - scaling parameter for beta RW (def = 2.93/sqrt(nvar))
//           s_alpha - scaling parameter for alpha RW (def = 2.93)
//           w - fractional weighting parameter (def = .1)
//           Vbeta0, Delta0 - initial guesses for parameters, if not supplied default values are used

  double ldiff, acc, unif, logalphac, oldlpostalpha, clpostalpha;
  int mkeep, rep;
  int nreg = regdata.size();
  int nz = Z.n_cols;
  int nvar = rootA.n_cols;  
  int nacceptbeta = 0;
  int nacceptalpha = 0; 
  
  mat Vbetainv = trans(rootA)*rootA;

  // allocate space for draws
  vec oldlpostbeta = zeros<vec>(nreg);
  vec clpostbeta = zeros<vec>(nreg);
  cube Betadraw = zeros<cube>(nreg, nvar, R/keep);
  vec alphadraw = zeros<vec>(R/keep);
  vec llike = zeros<vec>(R/keep);
  mat Vbetadraw = zeros<mat>(R/keep,nvar*nvar);
  mat Deltadraw = zeros<mat>(R/keep,nvar*nz);
  
  // convert regdata and hessdata Lists to std::vector of struct
  std::vector<moments> regdata_vector;
  moments regdatai_struct;
  List regdatai,hessi;

  // store vector with struct
  for (int reg = 0; reg<nreg; reg++){
    regdatai = regdata[reg];
    hessi = hessdata[reg];
  
    regdatai_struct.y = as<vec>(regdatai["y"]);
    regdatai_struct.X = as<mat>(regdatai["X"]);
    regdatai_struct.hess = as<mat>(hessi["hess"]);
    regdata_vector.push_back(regdatai_struct);    
  }

  if (nprint>0) startMcmcTimer();
  
  //  start main iteration loop
  for (rep = 0; rep < R; rep++){
    
    mat betabar = Z*Delta;
    
    // Draw betai
    for(int reg = 0; reg<nreg; reg++){
        vec betabari = trans(betabar(reg,span::all));
        mat betacvar = sbeta*solve(regdata_vector[reg].hess+Vbetainv,eye(nvar,nvar));
        mat betaroot = trans(chol(betacvar));
        vec betac = vectorise(Beta(reg,span::all)) + betaroot*vec(rnorm(nvar));
       
        oldlpostbeta[reg] = lpostbeta(alpha, trans(Beta(reg,span::all)), regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootA);
        clpostbeta[reg] = lpostbeta(alpha, betac, regdata_vector[reg].X, regdata_vector[reg].y, betabari, rootA);
        ldiff = clpostbeta[reg] - oldlpostbeta[reg];
        acc = exp(ldiff);
        if (acc > 1) acc = 1;    
        if(acc < 1) {unif=runif(1)[0];} else {unif=0;} //runif returns a NumericVector, so using [0] allows for conversion to double by extracting the first element
        if (unif <= acc){
          Beta(reg,span::all) = trans(betac);
          nacceptbeta = nacceptbeta + 1;
        }
    }
    
    // Draw alpha
    if (!fixalpha){
      logalphac = log(alpha) + alphacroot*rnorm(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
      oldlpostalpha = llnegbinpooled(regdata_vector,Beta,alpha)+(a-1)*log(alpha) - b*alpha;
      clpostalpha = llnegbinpooled(regdata_vector,Beta,exp(logalphac))+(a-1)*logalphac - b*exp(logalphac);
      ldiff = clpostalpha - oldlpostalpha;
      acc = exp(ldiff);
      if (acc > 1) acc = 1;    
      if(acc < 1) {unif=runif(1)[0];} else {unif=0;} //runif returns a NumericVector, so using [0] allows for conversion to double by extracting the first element
      if (unif <= acc){
        alpha = exp(logalphac);
        nacceptalpha = nacceptalpha + 1;
      }
    }  

    // Draw Vbeta and Delta using rmultireg
    List temp = rmultireg(Beta,Z,Deltabar,Adelta,nu,V);
    mat Vbeta = as<mat>(temp["Sigma"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    Vbetainv = solve(Vbeta,eye(nvar,nvar));
    rootA = chol(Vbetainv);
    Delta = as<mat>(temp["B"]);

    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);    
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      Betadraw.slice(mkeep-1) = Beta;
      alphadraw[mkeep-1] = alpha;
      Vbetadraw(mkeep-1,span::all) = trans(vectorise(Vbeta));
      Deltadraw(mkeep-1,span::all) = trans(vectorise(Delta));
      llike[mkeep-1] = llnegbinpooled(regdata_vector,Beta,alpha);
      } 
  }
  
  if (nprint>0) endMcmcTimer();
  
  return List::create(
    Named("llike") = llike,
    Named("Betadraw") = Betadraw,
    Named("alphadraw") = alphadraw,      
    Named("Vbetadraw") = Vbetadraw,
    Named("Deltadraw") = Deltadraw,
    Named("acceptrbeta") = nacceptbeta/(R*nreg*1.0)*100,
    Named("acceptralpha") = nacceptalpha/(R*1.0)*100);
}
