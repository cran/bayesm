#include "bayesm.h"
 
//FUNCTIONS SPECIFIC TO MAIN FUNCTION------------------------------------------------------
mat drawDelta(mat const& x,mat const& y,ivec const& z,std::vector<murooti> const& comps_vector,vec const& deltabar,mat const& Ad){

// Wayne Taylor 2/21/2015

// delta = vec(D)
//  given z and comps (z[i] gives component indicator for the ith observation, 
//   comps is a list of mu and rooti)
// y is n x p
// x is n x k
// y = xD' + U , rows of U are indep with covs Sigma_i given by z and comps

  int p = y.n_cols;
  int k = x.n_cols;
  int ncomp  = comps_vector.size();
  mat xtx = zeros<mat>(k*p,k*p);
  mat xty = zeros<mat>(p,k); //this is the unvecced version, reshaped after the sum
  
  //Create the index vectors, the colAll vectors are equal to span::all but with uvecs (as required by .submat)
  uvec colAlly(p), colAllx(k);
  for(int i = 0; i<p; i++) colAlly(i) = i;
  for(int i = 0; i<k; i++) colAllx(i) = i;
  
  //Loop through the components
  for(int compi = 0; compi<ncomp; compi++){
    
    //Create an index vector ind, to be used like y[ind,]
    uvec ind = find(z == (compi+1));
  
    //If there are observations in this component
    if(ind.size()>0){
      mat yi = y.submat(ind,colAlly);
      mat xi = x.submat(ind,colAllx);
      
      murooti compsi_struct = comps_vector[compi];
      yi.each_row() -= trans(compsi_struct.mu); //the subtraction operation is repeated on each row of yi
      mat sigi = compsi_struct.rooti*trans(compsi_struct.rooti);
      xtx = xtx + kron(trans(xi)*xi,sigi);
      xty = xty + (sigi * (trans(yi)*xi));
    }
  }
  xty.reshape(xty.n_rows*xty.n_cols,1);
  
  //vec(t(D)) ~ N(V^{-1}(xty + Ad*deltabar),V^{-1}) where V = (xtx+Ad)
  // compute the inverse of xtx+Ad
  mat ucholinv = solve(trimatu(chol(xtx+Ad)), eye(k*p,k*p)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  mat Vinv = ucholinv*trans(ucholinv);
  
  return(Vinv*(xty+Ad*deltabar) + trans(chol(Vinv))*as<vec>(rnorm(deltabar.size())));
}

DPOut rDPGibbs1(mat y, lambda lambda_struct, std::vector<murooti> thetaStar_vector, int maxuniq, ivec indic, 
            vec q0v, double alpha, priorAlpha const& priorAlpha_struct, int gridsize, List const& lambda_hyper){

// Wayne Taylor 2/21/2015

//revision history:
//created from rDPGibbs by Rossi 3/08

//do one draw of DP Gibbs sampler with normal base

//Model:
//  y_i ~ N(y|thetai)
//  thetai|G ~ G
//  G|lambda,alpha ~ DP(G|G0(lambda),alpha)

//Priors:
//  alpha: starting value
//  lambda:
//    G0 ~ N(mubar,Sigma (x) Amu^-1)
//    mubar=vec(mubar)
//    Sigma ~ IW(nu,nu*V) V=v*I  note: mode(Sigma)=nu/(nu+2)*v*I
//    mubar=0
//    amu is uniform on grid specified by alim
//    nu is log uniform, nu=d-1+exp(Z) z is uniform on seq defined bvy nulim
//    v is uniform on sequence specificd by vlim

//  priorAlpha_struct:
//    alpha ~ (1-(alpha-alphamin)/(alphamax-alphamin))^power
//    alphamin=exp(digamma(Istarmin)-log(gamma+log(N)))
//    alphamax=exp(digamma(Istarmax)-log(gamma+log(N)))
//    gamma= .5772156649015328606

//output:
//  ind - vector of indicators for which observations are associated with which comp in thetaStar
//  thetaStar - list of unique normal component parms
//  lambda  - list of of (a,nu,V)
//  alpha 
//  thetaNp1 - one draw from predictive given thetaStar, lambda,alphama

  int n = y.n_rows;
  int dimy = y.n_cols;
  int nunique, indsize, indp, probssize;
  vec probs;
  uvec ind;
  mat ydenmat;
  uvec spanall(dimy); for(int i = 0; i<dimy ; i++) spanall[i] = i; //creates a uvec of [0,1,...,dimy-1]
  thetaStarIndex thetaStarDrawOut_struct;
  std::vector<murooti> new_utheta_vector(1), thetaNp1_vector(1);
  murooti thetaNp10_struct, outGD_struct;

  for(int rep = 0; rep<1; rep++) { //note we only do one loop!
    
    q0v = q0(y,lambda_struct);
   
    nunique = thetaStar_vector.size();
  
    if(nunique > maxuniq) stop("maximum number of unique thetas exceeded");
   
    //ydenmat is a length(thetaStar) x n array of density values given f(y[j,] | thetaStar[[i]]
    //  note: due to remix step (below) we must recompute ydenmat each time!
    ydenmat = zeros<mat>(maxuniq,n);
                     
    ydenmat(span(0,nunique-1),span::all) = yden(thetaStar_vector,y);
  
    thetaStarDrawOut_struct = thetaStarDraw(indic, thetaStar_vector, y, ydenmat, q0v, alpha, lambda_struct, maxuniq);
    thetaStar_vector = thetaStarDrawOut_struct.thetaStar_vector;
    indic = thetaStarDrawOut_struct.indic;
    nunique = thetaStar_vector.size();
  
    //thetaNp1 and remix
    probs = zeros<vec>(nunique+1);
    for(int j = 0; j < nunique; j++){
      ind = find(indic == (j+1));
      indsize = ind.size();
      probs[j] = indsize/(alpha + n + 0.0);
      new_utheta_vector[0] = thetaD(y(ind,spanall),lambda_struct);
      thetaStar_vector[j] = new_utheta_vector[0];
    }
                  
    probs[nunique] = alpha/(alpha+n+0.0);
    indp = rmultinomF(probs);
    probssize = probs.size();
    if(indp == probssize) {
      outGD_struct = GD(lambda_struct);
      thetaNp10_struct.mu = outGD_struct.mu;
      thetaNp10_struct.rooti = outGD_struct.rooti;
      thetaNp1_vector[0] = thetaNp10_struct;
    } else {
      outGD_struct = thetaStar_vector[indp-1];
      thetaNp10_struct.mu = outGD_struct.mu;
      thetaNp10_struct.rooti = outGD_struct.rooti;
      thetaNp1_vector[0] = thetaNp10_struct;
    }
  
    //draw alpha
    alpha = alphaD(priorAlpha_struct,nunique,gridsize);
  
    //draw lambda
    lambda_struct = lambdaD(lambda_struct,thetaStar_vector,lambda_hyper["alim"],lambda_hyper["nulim"],lambda_hyper["vlim"],gridsize);
  }

  //note indic is the vector of indicators for each obs correspond to which thetaStar
  DPOut out_struct;
    out_struct.thetaStar_vector = thetaStar_vector;
    out_struct.thetaNp1_vector = thetaNp1_vector;
    out_struct.alpha = alpha;
    out_struct.lambda_struct = lambda_struct;
    out_struct.indic = indic;

  return(out_struct);
}

//MAIN FUNCTION-------------------------------------------------------------------------------------
//[[Rcpp::export]]
List rhierMnlDP_rcpp_loop(int R, int keep, int nprint,
                          List const& lgtdata, mat const& Z,
                          vec const& deltabar, mat const& Ad, List const& PrioralphaList, List const& lambda_hyper,
                          bool drawdelta, int nvar, mat oldbetas, double s,
                          int maxuniq, int gridsize,
                          double BayesmConstantA, int BayesmConstantnuInc, double BayesmConstantDPalpha){

// Wayne Taylor 2/21/2015

  //Initialize variable placeholders
  int mkeep, Istar;
  vec betabar, q0v;
  mat rootpi, ucholinv, incroot, V;
  List compdraw(R/keep), nmix;
  DPOut mgout_struct;
  mnlMetropOnceOut metropout_struct;
  murooti thetaStarLgt_struct;
  
  int nz = Z.n_cols;
  int nlgt = lgtdata.size();
  
  // convert List to std::vector of struct
  List lgtdatai;
  std::vector<moments> lgtdata_vector;
  moments lgtdatai_struct;
  for (int lgt = 0; lgt<nlgt; lgt++){
    lgtdatai = lgtdata[lgt];
    
    lgtdatai_struct.y = as<vec>(lgtdatai["y"]);
    lgtdatai_struct.X = as<mat>(lgtdatai["X"]);
    lgtdatai_struct.hess = as<mat>(lgtdatai["hess"]);
    lgtdata_vector.push_back(lgtdatai_struct);    
  }
  
  //initialize indicator vector, delta, thetaStar, thetaNp10, alpha, oldprob
  ivec indic = ones<ivec>(nlgt);
  
  mat olddelta;
  if (drawdelta) olddelta = zeros<vec>(nz*nvar);
   
  std::vector<murooti> thetaStar_vector(1);
  murooti thetaNp10_struct, thetaStar0_struct;
    thetaStar0_struct.mu = zeros<vec>(nvar);
    thetaStar0_struct.rooti = eye(nvar,nvar);
  thetaStar_vector[0] = thetaStar0_struct;
  
  double alpha = BayesmConstantDPalpha;

  //fix oldprob (only one comp)
  double oldprob = 1.0;

  //convert Prioralpha from List to struct
  priorAlpha priorAlpha_struct;
    priorAlpha_struct.power = PrioralphaList["power"];
    priorAlpha_struct.alphamin = PrioralphaList["alphamin"];
    priorAlpha_struct.alphamax = PrioralphaList["alphamax"];
    priorAlpha_struct.n = PrioralphaList["n"];

 //initialize lambda
  lambda lambda_struct;
    lambda_struct.mubar = zeros<vec>(nvar);
    lambda_struct.Amu = BayesmConstantA;
    lambda_struct.nu = nvar+BayesmConstantnuInc;
    lambda_struct.V = lambda_struct.nu*eye(nvar,nvar);
  
  //allocate space for draws
  mat Deltadraw(1,1); if(drawdelta) Deltadraw.zeros(R/keep, nz*nvar);//enlarge Deltadraw only if the space is required
  cube betadraw(nlgt, nvar, R/keep);
  vec probdraw = zeros<vec>(R/keep);
  vec oldll = zeros<vec>(nlgt);
  vec loglike = zeros<vec>(R/keep);
  vec Istardraw = zeros<vec>(R/keep);
  vec alphadraw = zeros<vec>(R/keep);
  vec nudraw = zeros<vec>(R/keep);
  vec vdraw = zeros<vec>(R/keep);
  vec adraw = zeros<vec>(R/keep);
  
  if (nprint>0) startMcmcTimer();
  
  //start main iteration loop
  for(int rep = 0; rep<R; rep++) {
    
    //first draw comps,indic,p | {beta_i}, delta
    //  indic,p need initialization comps is drawn first in sub-Gibbs
    if(drawdelta){
      olddelta.reshape(nvar,nz);
      mgout_struct = rDPGibbs1(oldbetas-Z*trans(olddelta),lambda_struct,thetaStar_vector,maxuniq,indic,q0v,alpha,priorAlpha_struct,gridsize,lambda_hyper);
    } else {
      mgout_struct = rDPGibbs1(oldbetas,lambda_struct,thetaStar_vector,maxuniq,indic,q0v,alpha,priorAlpha_struct,gridsize,lambda_hyper);
    }
  
    indic = mgout_struct.indic;
    lambda_struct = mgout_struct.lambda_struct;
    alpha = mgout_struct.alpha;
    thetaStar_vector = mgout_struct.thetaStar_vector;
    Istar = thetaStar_vector.size();

    //now draw delta | {beta_i}, ind, comps
    if(drawdelta) {olddelta = drawDelta(Z,oldbetas,indic,thetaStar_vector,deltabar,Ad);}
  
    //loop over all lgt equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
    for (int lgt=0; lgt<nlgt; lgt++){
      thetaStarLgt_struct = thetaStar_vector[indic[lgt]-1];
      rootpi = thetaStarLgt_struct.rooti;
      
      //note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
      if(drawdelta){
        olddelta.reshape(nvar,nz);
        betabar = thetaStarLgt_struct.mu + olddelta * trans(Z(lgt,span::all));
      } else {
        betabar = thetaStarLgt_struct.mu;
      }
      
      if (rep == 0) oldll[lgt] = llmnl(vectorise(oldbetas(lgt,span::all)),lgtdata_vector[lgt].y,lgtdata_vector[lgt].X);
      
      //compute inc.root
      ucholinv = solve(trimatu(chol(lgtdata_vector[lgt].hess+rootpi*trans(rootpi))), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
      incroot = chol(ucholinv*trans(ucholinv));
      
      metropout_struct = mnlMetropOnce(lgtdata_vector[lgt].y,lgtdata_vector[lgt].X,vectorise(oldbetas(lgt,span::all)),
                                           oldll[lgt],s,incroot,betabar,rootpi);
      
      oldbetas(lgt,span::all) = trans(metropout_struct.betadraw);
      oldll[lgt] = metropout_struct.oldll;   
    }
  
    //print time to completion and draw # every nprint'th draw
    if (nprint>0) if((rep+1)%nprint==0) infoMcmcTimer(rep, R);
      
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw.slice(mkeep-1) = oldbetas;
      probdraw[mkeep-1] = oldprob;
      alphadraw[mkeep-1] = alpha;
      Istardraw[mkeep-1] = Istar;
      adraw[mkeep-1] = lambda_struct.Amu;
      nudraw[mkeep-1] = lambda_struct.nu;
      V = lambda_struct.V;
      vdraw[mkeep-1] = V(0,0)/(lambda_struct.nu+0.0);
      loglike[mkeep-1] = sum(oldll);
      if(drawdelta) Deltadraw(mkeep-1, span::all) = trans(vectorise(olddelta));
      thetaNp10_struct = mgout_struct.thetaNp1_vector[0];
      //we have to convert to a NumericVector for the plotting functions to work
      compdraw[mkeep-1] = List::create(List::create(Named("mu") = NumericVector(thetaNp10_struct.mu.begin(),thetaNp10_struct.mu.end()),Named("rooti") = thetaNp10_struct.rooti));
    }
  }
  
  if (nprint>0) endMcmcTimer();
  
  nmix = List::create(Named("probdraw") = probdraw,
    			  Named("zdraw") = R_NilValue, //sets the value to NULL in R
					  Named("compdraw") = compdraw);
	
  if(drawdelta){
    return(List::create(
      Named("Deltadraw") = Deltadraw,
      Named("betadraw") = betadraw,
      Named("nmix") = nmix,
      Named("alphadraw") = alphadraw,
      Named("Istardraw") = Istardraw,
      Named("adraw") = adraw,
      Named("nudraw") = nudraw,
      Named("vdraw") = vdraw,
      Named("loglike") = loglike));
	} else {
    return(List::create(
      Named("betadraw") = betadraw,
      Named("nmix") = nmix,
      Named("alphadraw") = alphadraw,
      Named("Istardraw") = Istardraw,
      Named("adraw") = adraw,
      Named("nudraw") = nudraw,
      Named("vdraw") = vdraw,
      Named("loglike") = loglike));
	}
}
