#include "bayesm.h"
 
//[[Rcpp::export]]
List rDPGibbs_rcpp_loop(int R, int keep, int nprint,
                        mat y, List const& lambda_hyper, bool SCALE, int maxuniq, List const& PrioralphaList, int gridsize,
                        double BayesmConstantA, int BayesmConstantnuInc, double BayesmConstantDPalpha) {

// Wayne Taylor 2/4/2015

  int dimy = y.n_cols;
  int n = y.n_rows;
  
  //initialize indic, thetaStar, thetaNp1
  ivec indic = ones<ivec>(n);
  
  std::vector<murooti> thetaStar_vector(1);
  murooti thetaStar0_struct;
    thetaStar0_struct.mu = zeros<vec>(dimy);
    thetaStar0_struct.rooti = eye(dimy,dimy);
  thetaStar_vector[0] = thetaStar0_struct;

  //convert Prioralpha from List to struct
  priorAlpha priorAlpha_struct;
    priorAlpha_struct.power = PrioralphaList["power"];
    priorAlpha_struct.alphamin = PrioralphaList["alphamin"];
    priorAlpha_struct.alphamax = PrioralphaList["alphamax"];
    priorAlpha_struct.n = PrioralphaList["n"];

 //initialize lambda
  lambda lambda_struct;
    lambda_struct.mubar = zeros<vec>(dimy);
    lambda_struct.Amu = BayesmConstantA;
    lambda_struct.nu = dimy+BayesmConstantnuInc;
    lambda_struct.V = lambda_struct.nu*eye(dimy,dimy);

  //initialize alpha
  double alpha = BayesmConstantDPalpha;
  
  //intialize remaining variables
  thetaStarIndex thetaStarDrawOut_struct;
  std::vector<murooti> new_utheta_vector(1), thetaNp1_vector(1);
  murooti thetaNp10_struct, out_struct;
  mat ydenmat;
  vec q0v, probs;
  uvec ind;
  int nunique, indsize;
  uvec spanall(dimy); for(int i = 0; i<dimy ; i++) spanall[i] = i; //creates a uvec of [0,1,...,dimy-1]
  double nu;

  //allocate storage
  vec alphadraw = zeros<vec>(R/keep);
  vec Istardraw = zeros<vec>(R/keep);
  vec adraw = zeros<vec>(R/keep);
  vec nudraw = zeros<vec>(R/keep);
  vec vdraw = zeros<vec>(R/keep);
  List thetaNp1draw(R/keep);
  imat inddraw = zeros<imat>(R/keep,n);

  //do scaling
  rowvec dvec, ybar;
  if(SCALE){
    dvec = 1/sqrt(var(y,0,0)); //norm_type=0 performs normalisation using N-1, dim=0 is by column
    ybar = mean(y,0);
    y.each_row() -= ybar; //subtract ybar from each row
    y.each_row() %= dvec; //divide each row by dvec
  } 
  //note on scaling
  //we model scaled y, z_i=D(y_i-ybar)   D=diag(1/sigma1, ..., 1/sigma_dimy)
  
  //if p_z= 1/R sum(phi(z|mu,Sigma))
  // p_y=1/R sum(phi(y|D^-1mu+ybar,D^-1SigmaD^-1)
  // rooti_y=Drooti_z
  
  //you might want to use quantiles instead, like median and (10,90)

  // start main iteration loop
  int mkeep = 0;
  
  if(nprint>0) startMcmcTimer();

  for(int rep = 0; rep<R; rep++) {
    
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
    int ind = rmultinomF(probs);
    int probssize = probs.size();
    if(ind == probssize) {
      out_struct = GD(lambda_struct);
      thetaNp10_struct.mu = out_struct.mu;
      thetaNp10_struct.rooti = out_struct.rooti;
      thetaNp1_vector[0] = thetaNp10_struct;
    } else {
      out_struct = thetaStar_vector[ind-1];
      thetaNp10_struct.mu = out_struct.mu;
      thetaNp10_struct.rooti = out_struct.rooti;
      thetaNp1_vector[0] = thetaNp10_struct;
    }
  
    //draw alpha
    alpha = alphaD(priorAlpha_struct,nunique,gridsize);
  
    //draw lambda
    lambda_struct = lambdaD(lambda_struct,thetaStar_vector,lambda_hyper["alim"],lambda_hyper["nulim"],lambda_hyper["vlim"],gridsize);
  
    //print time to completion
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
   
    //save every keepth draw
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      alphadraw[mkeep-1] = alpha;
      Istardraw[mkeep-1] = nunique;
      adraw[mkeep-1] = lambda_struct.Amu;
      nu = lambda_struct.nu;
      nudraw[mkeep-1] = nu;
      mat V = lambda_struct.V;
      vdraw[mkeep-1] = V(0,0)/(nu+0.0);
      inddraw(mkeep-1,span::all) = trans(indic);
      
      thetaNp10_struct = thetaNp1_vector[0];
      if(SCALE){
        thetaNp10_struct.mu = thetaNp10_struct.mu/trans(dvec)+trans(ybar);
        thetaNp10_struct.rooti = diagmat(dvec)*thetaNp10_struct.rooti;
      }

      //here we put the draws into the list of lists of list format useful for finite mixture of normals utilities
      //we have to convetr to a NumericVector for the plotting functions to work
      thetaNp1draw[mkeep-1] = List::create(List::create(Named("mu") = NumericVector(thetaNp10_struct.mu.begin(),thetaNp10_struct.mu.end()),Named("rooti") = thetaNp10_struct.rooti));
     }  
  }

  if(nprint>0) endMcmcTimer();
  
  return List::create(
    Named("inddraw") = inddraw,
    Named("thetaNp1draw") = thetaNp1draw,
    Named("alphadraw") = alphadraw,
    Named("Istardraw") = Istardraw,
    Named("adraw") = adraw,
    Named("nudraw") = nudraw,
    Named("vdraw") = vdraw);
}
