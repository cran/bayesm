#include "bayesm.h"

//FUNCTION SPECIFIC TO MAIN FUNCTION------------------------------------------------------
//[[Rcpp::export]]
double llmnl_con(vec const& betastar, vec const& y, mat const& X, vec const& SignRes = NumericVector::create(0)){
  
  // Wayne Taylor 7/8/2016
  
  // Evaluates log-likelihood for the multinomial logit model WITH SIGN CONSTRAINTS
  // NOTE: this is exported only because it is used in the shell .R function, it will not be available to users
  
  //Reparameterize betastar to beta to allow for sign restrictions
  vec beta = betastar;
  
  //The default SignRes vector is a single element vector containing a zero
  //any() returns true if any elements of SignRes are non-zero
  if(any(SignRes)){ 
    uvec signInd = find(SignRes != 0);
    beta.elem(signInd) = SignRes.elem(signInd) % exp(beta.elem(signInd));  //% performs element-wise multiplication
  }
  
  int n = y.size();
  int j = X.n_rows/n;
  mat Xbeta = X*beta;
  
  vec xby = zeros<vec>(n);
  vec denom = zeros<vec>(n);
  
  for(int i = 0; i<n;i++){      
    for(int p=0;p<j;p++) denom[i]=denom[i]+exp(Xbeta[i*j+p]);
    xby[i] = Xbeta[i*j+y[i]-1];
  }
  
  return(sum(xby - log(denom)));
}

//mnlRwMetropOnce=
//function(y,X,oldbeta,oldll,s,inc.root,betabar,rootpi){ 
//#
//# function to execute rw metropolis for the MNL
//# y is n vector with element = 1,...,j indicating which alt chosen
//# X is nj x k matrix of xvalues for each of j alt on each of n occasions
//# RW increments are N(0,s^2*t(inc.root)%*%inc.root)
//# prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
//#  inc.root, rootpi are upper triangular
//#  this means that we are using the UL decomp of Sigma^-1 for prior 
//# oldbeta is the current
//     stay=0
//     betac=oldbeta + s*t(inc.root)%*%(matrix(rnorm(ncol(X)),ncol=1))
//     cll=llmnl(betac,y,X)
//     clpost=cll+lndMvn(betac,betabar,rootpi)
//     ldiff=clpost-oldll-lndMvn(oldbeta,betabar,rootpi)
//     alpha=min(1,exp(ldiff))
//     if(alpha < 1) {unif=runif(1)} else {unif=0}
//     if (unif <= alpha)
//             {betadraw=betac; oldll=cll}
//           else
//             {betadraw=oldbeta; stay=1}
//return(list(betadraw=betadraw,stay=stay,oldll=oldll))
//}

mnlMetropOnceOut mnlMetropOnce_con(vec const& y, mat const& X, vec const& oldbeta, 
                                            double oldll,double s, mat const& incroot, 
                                            vec const& betabar, mat const& rootpi,vec const& SignRes = NumericVector::create(2)){ 
  // Wayne Taylor 10/01/2014
  
  // function to execute rw metropolis for the MNL
  // y is n vector with element = 1,...,j indicating which alt chosen
  // X is nj x k matrix of xvalues for each of j alt on each of n occasions
  // RW increments are N(0,s^2*t(inc.root)%*%inc.root)
  // prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
  //  inc.root, rootpi are upper triangular
  //  this means that we are using the UL decomp of Sigma^-1 for prior 
  // oldbeta is the current
  
  
  mnlMetropOnceOut out_struct;
  
  double unif;
  vec betadraw, alphaminv;
  
  int stay = 0;
  vec betac = oldbeta + s*trans(incroot)*as<vec>(rnorm(X.n_cols));
  double cll = llmnl_con(betac,y,X,SignRes);
  double clpost = cll+lndMvn(betac,betabar,rootpi);
  double ldiff = clpost-oldll-lndMvn(oldbeta,betabar,rootpi);
  alphaminv << 1 << exp(ldiff);
  double alpha = min(alphaminv);
  
  if(alpha < 1) {
    unif = as_scalar(vec(runif(1)));
  } else { 
    unif=0;}
  if (unif <= alpha) {
    betadraw = betac;
    oldll = cll;
  } else {
    betadraw = oldbeta;
    stay = 1;
  }
  
  out_struct.betadraw = betadraw;
  out_struct.stay = stay;  
  out_struct.oldll = oldll;
  
  return (out_struct);
}

//MAIN FUNCTION-------------------------------------------------------------------------------------

//[[Rcpp::export]]
List rhierMnlRwMixture_rcpp_loop(List const& lgtdata, mat const& Z,
                                  vec const& deltabar, mat const& Ad, mat const& mubar, mat const& Amu,
                                  double nu, mat const& V, double s,
                                  int R, int keep, int nprint, bool drawdelta,
                                  mat olddelta,  vec const& a, vec oldprob, mat oldbetas, vec ind, vec const& SignRes){

// Wayne Taylor 10/01/2014

  int nlgt = lgtdata.size();
  int nvar = V.n_cols;
  int nz = Z.n_cols;
  
  mat rootpi, betabar, ucholinv, incroot;
  int mkeep;
  mnlMetropOnceOut metropout_struct;
  List lgtdatai, nmix;
  
  // convert List to std::vector of struct
  std::vector<moments> lgtdata_vector;
  moments lgtdatai_struct;
  for (int lgt = 0; lgt<nlgt; lgt++){
    lgtdatai = lgtdata[lgt];
    
    lgtdatai_struct.y = as<vec>(lgtdatai["y"]);
    lgtdatai_struct.X = as<mat>(lgtdatai["X"]);
    lgtdatai_struct.hess = as<mat>(lgtdatai["hess"]);
    lgtdata_vector.push_back(lgtdatai_struct);    
  }
    
  // allocate space for draws
  vec oldll = zeros<vec>(nlgt);
  cube betadraw(nlgt, nvar, R/keep);
  mat probdraw(R/keep, oldprob.size());
  vec loglike(R/keep);
  mat Deltadraw(1,1); if(drawdelta) Deltadraw.zeros(R/keep, nz*nvar);//enlarge Deltadraw only if the space is required
  List compdraw(R/keep);
  
  if (nprint>0) startMcmcTimer();
    
  for (int rep = 0; rep<R; rep++){
    
    //first draw comps,ind,p | {beta_i}, delta
    // ind,p need initialization comps is drawn first in sub-Gibbs
    List mgout;
    if(drawdelta) {
      olddelta.reshape(nvar,nz);
      mgout = rmixGibbs (oldbetas-Z*trans(olddelta),mubar,Amu,nu,V,a,oldprob,ind);
    } else {
      mgout = rmixGibbs(oldbetas,mubar,Amu,nu,V,a,oldprob,ind);
    }
    
    List oldcomp = mgout["comps"];
    oldprob = as<vec>(mgout["p"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    ind = as<vec>(mgout["z"]);
    
    //now draw delta | {beta_i}, ind, comps
    if(drawdelta) olddelta = drawDelta(Z,oldbetas,ind,oldcomp,deltabar,Ad);
    
    //loop over all LGT equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
      for(int lgt = 0; lgt<nlgt; lgt++){
        List oldcomplgt = oldcomp[ind[lgt]-1];
        rootpi = as<mat>(oldcomplgt[1]);
        
        //note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
        if(drawdelta){
          olddelta.reshape(nvar,nz);
          betabar = as<vec>(oldcomplgt[0])+olddelta*vectorise(Z(lgt,span::all));
        } else {
          betabar = as<vec>(oldcomplgt[0]);
        }
        
        if (rep == 0) oldll[lgt] = llmnl_con(vectorise(oldbetas(lgt,span::all)),lgtdata_vector[lgt].y,lgtdata_vector[lgt].X,SignRes);
        
        //compute inc.root
        ucholinv = solve(trimatu(chol(lgtdata_vector[lgt].hess+rootpi*trans(rootpi))), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
        incroot = chol(ucholinv*trans(ucholinv));
                
        metropout_struct = mnlMetropOnce_con(lgtdata_vector[lgt].y,lgtdata_vector[lgt].X,vectorise(oldbetas(lgt,span::all)),
                                         oldll[lgt],s,incroot,betabar,rootpi,SignRes);
         
         oldbetas(lgt,span::all) = trans(metropout_struct.betadraw);
         oldll[lgt] = metropout_struct.oldll;  
      }
      
    //print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw.slice(mkeep-1) = oldbetas;
      probdraw(mkeep-1, span::all) = trans(oldprob);
      loglike[mkeep-1] = sum(oldll);
      if(drawdelta) Deltadraw(mkeep-1, span::all) = trans(vectorise(olddelta));
      compdraw[mkeep-1] = oldcomp;
    }
  }
  
  if (nprint>0) endMcmcTimer();
  
  nmix = List::create(Named("probdraw") = probdraw,
    		  Named("zdraw") = R_NilValue, //sets the value to NULL in R
				  Named("compdraw") = compdraw);
  
  //ADDED FOR CONSTRAINTS
  //If there are sign constraints, return f(betadraws) as "betadraws"
  //conStatus will be set to true if SignRes has any non-zero elements
  bool conStatus = any(SignRes);
  
  if(conStatus){
    int SignResSize = SignRes.size();
    
    //loop through each sign constraint
    for(int i = 0;i < SignResSize; i++){
      
      //if there is a constraint loop through each slice of betadraw
      if(SignRes[i] != 0){
        for(int s = 0;s < R/keep; s++){
          betadraw(span(),span(i),span(s)) = SignRes[i] * exp(betadraw(span(),span(i),span(s)));
        }
      }
      
    }//end loop through SignRes
  }
  
  if(drawdelta){
    return(List::create(
        Named("Deltadraw") = Deltadraw,
        Named("betadraw") = betadraw,
        Named("nmix") = nmix,
        Named("loglike") = loglike,
        Named("SignRes") = SignRes));  
  } else {
    return(List::create(
        Named("betadraw") = betadraw,
        Named("nmix") = nmix,
        Named("loglike") = loglike,
        Named("SignRes") = SignRes));
  }
  
}
