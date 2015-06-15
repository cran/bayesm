#include "bayesm.h"
 
//FUNCTIONS SPECIFIC TO MAIN FUNCTION------------------------------------------------------
struct ytxtxtd{
  vec yt;
  vec xt;
  mat xtd;
};

ytxtxtd get_ytxt(vec const& y, mat const& z, mat const& delta, mat const& x, mat const& w,
                              int ncomp,ivec const& indic, std::vector<murooti> const& thetaStar_vector){

  // Wayne Taylor 3/14/2015
  
  int dimz = z.n_cols;
  int dimx = x.n_cols;
  
  //variable type initializaion
  double sig;
  mat wk, zk, xk, rooti, Sigma, xt;
  vec yk, mu, e1, ee2, yt;
  uvec ind, colAllw, colAllz(dimz), colAllx(dimx);

  //Create the index vectors, the colAll vectors are equal to span::all but with uvecs (as required by .submat)
  for(int i = 0; i<dimz; i++) colAllz(i) = i;
  for(int i = 0; i<dimx; i++) colAllx(i) = i;

  bool isw = false;
  if(!w.is_empty()){
    isw = true;
    int ncolw = w.n_cols;
    uvec colAllw(ncolw);
    for(int i = 0; i<ncolw; i++) colAllw(i) = i;
  }
  
  for (int k = 0; k < ncomp; k++){ 
    
    //Create an index vector ind, to be used like y[ind,]
    ind = find(indic == (k+1));
  
    //If there are observations in this component
    if(ind.size()>0){
      
      if(isw) wk = w.submat(ind,colAllw);
      zk = z.submat(ind,colAllz);
      yk = y(ind);
      xk = x.submat(ind,colAllx);
      
      murooti thetaStark_struct = thetaStar_vector[k];
      mu = thetaStark_struct.mu;
      rooti = thetaStark_struct.rooti;
      
      Sigma = solve(rooti,eye(2,2));
      Sigma = trans(Sigma)*Sigma;
    
      e1 = xk-zk*delta;  
      ee2 = mu[1] + (Sigma(0,1)/Sigma(0,0))*(e1-mu[0]);
      sig = sqrt(Sigma(1,1)-pow(Sigma(0,1),2.0)/Sigma(0,0));
      yt = join_cols(yt,(yk-ee2)/sig); //analogous to rbind()
      
      if(isw) {
        xt = join_cols(xt,join_rows(xk,wk)/sig);
      } else {
        xt = join_cols(xt,xk/sig);
      }
    }
  }
  
  ytxtxtd out_struct;
    out_struct.yt = yt;
    out_struct.xt = xt;
    
  return(out_struct);
}

ytxtxtd get_ytxtd(vec const& y, mat const& z, double beta, vec const& gamma, mat const& x, mat const& w,
                              int ncomp, ivec const& indic,std::vector<murooti> const& thetaStar_vector, int dimd){

  // Wayne Taylor 3/14/2015

  int dimx = x.n_cols;

  //variable type initializaion
  int indsize, indicsize;
  vec zveck, yk, mu, ytk, u, yt;
  mat C, wk, zk, xk, rooti, Sigma, B, L, Li, z2, zt1, zt2, xtd;
  uvec colAllw, colAllz(dimd), colAllx(dimx), ind, seqindk, negseqindk;

  //Create index vectors (uvec) for submatrix views
  indicsize = indic.size();
  //here the uvecs are declared once, and within each loop the correctly sized vector is extracted as needed
  uvec seqind(indicsize);for(int i = 0;i<indicsize;i++){seqind[i] = i*2;} //element 0,2,4,. . .
  uvec negseqind(indicsize);for(int i = 0;i<indicsize;i++){negseqind[i] = (i*2)+1;} //element 1,3,5,...
  
  //colAll vectors are equal to span::all but with uvecs (as required by .submat)
  for(int i = 0; i<dimd; i++) colAllz(i) = i;
  for(int i = 0; i<dimx; i++) colAllx(i) = i;
  
  bool isw = false;
  if(!w.is_empty()){
    isw = true;
    int ncolw = w.n_cols;
    uvec colAllw(ncolw);
    for(int i = 0; i<ncolw; i++) colAllw(i) = i;
  }
  
  C = eye(2,2); C(1,0) = beta;
  
  for(int k = 0;k<ncomp;k++){
      //Create an index vector ind, to be used like y[ind,]
    ind = find(indic == (k+1));
    indsize = ind.size();
  
    //If there are observations in this component
    if(indsize>0){
  
      mat xtdk(2*indsize,dimd);    
      
      //extract the properly sized vector section
      seqindk = seqind.subvec(0,indsize-1);
      negseqindk = negseqind.subvec(0,indsize-1);
      
      if(isw) wk = w.submat(ind,colAllw);
      zk = z.submat(ind,colAllz);
      zveck = vectorise(trans(zk));
      yk = y(ind);
      xk = x.submat(ind,colAllx);
  
      murooti thetaStark_struct = thetaStar_vector[k];
      mu = thetaStark_struct.mu;
      rooti = thetaStark_struct.rooti;

      Sigma = solve(rooti,eye(2,2));
      Sigma = trans(Sigma)*Sigma;
      
      B = C*Sigma*trans(C);
      L = trans(chol(B));
      Li = solve(trimatl(L),eye(2,2)); // L is lower triangular, trimatl interprets the matrix as lower triangular and makes solve more efficient
      if(isw) {
        u = vectorise(yk-wk*gamma-mu[1]-beta*mu[0]);
      } else {
        u = vectorise(yk-mu[1]-beta*mu[0]);
      }
      
      ytk = vectorise(Li * join_cols(trans(xk-mu[0]),trans(u)));
      
      z2 = trans(join_rows(zveck,beta*zveck)); //join_rows is analogous to cbind()
      z2 = Li*z2;
      zt1 = z2(0,span::all);
      zt2 = z2(1,span::all);
      
      zt1.reshape(dimd,indsize);
      zt1 = trans(zt1);
      zt2.reshape(dimd,indsize);
      zt2=trans(zt2);
      
      xtdk(seqindk,colAllz) = zt1;
      xtdk(negseqindk,colAllz) = zt2;
      
      yt = join_cols(yt,ytk);
      xtd = join_cols(xtd,xtdk);    
    }
  }

  ytxtxtd out_struct;
    out_struct.yt = yt;
    out_struct.xtd = xtd;
    
  return(out_struct);  
}

DPOut rthetaDP(int maxuniq, double alpha, lambda lambda_struct, priorAlpha const& priorAlpha_struct, 
                              std::vector<murooti> thetaStar_vector, ivec indic, vec const& q0v, mat const& y, int gridsize,
                              List lambda_hyper){
 
  // Wayne Taylor 3/14/2015

//  function to make one draw from DP process 

//  P. Rossi 1/06
//  added draw of alpha 2/06
//  removed lambdaD,etaD and function arguments 5/06
//  removed thetaStar argument to .Call and creation of newthetaStar 7/06
//  removed q0 computations as eta is not drawn  7/06
//  changed for new version of thetadraw and removed calculation of thetaStar before
//    .Call  7/07

//      y(i) ~ f(y|theta[[i]],eta)
//      theta ~ DP(alpha,G(lambda))

//output:
//   list with components:
//      thetaDraws: list, [[i]] is a list of the ith draw of the n theta's
//                  where n is the length of the input theta and nrow(y)
//      thetaNp1Draws: list, [[i]] is ith draw of theta_{n+1}
//args:
//   maxuniq: the maximum number of unique thetaStar values -- an error will be raised
//            if this is exceeded
//   alpha,lambda: starting values (or fixed DP prior values if not drawn).
//   Prioralpha: list of hyperparms of alpha prior
//   theta: list of starting value for theta's
//   thetaStar: list of unique values of theta, thetaStar[[i]]
//   indic:  n vector of indicator for which unique theta (in thetaStar)
//   y: is a matrix nxk
//         thetaStar: list of unique values of theta, thetaStar[[i]]
//   q0v:a double vector with the same number of rows as y, giving \Int f(y(i)|theta,eta) dG_{lambda}(theta).

  int n = y.n_rows;
  int dimy = y.n_cols;
  
  //variable type initializaion
  int nunique, indsize, indp, probssize;
  vec probs;
  uvec ind;
  mat ydenmat;
  uvec spanall(dimy); for(int i = 0; i<dimy ; i++) spanall[i] = i; //creates a uvec of [0,1,...,dimy-1]
  thetaStarIndex thetaStarDrawOut_struct;
  std::vector<murooti> new_utheta(1), thetaNp1_vector(1);
  murooti thetaNp10_struct, outGD;
  
  vec p(n);
  p[n-1] =  alpha/(alpha+(n-1));
  for(int i = 0; i<(n-1); i++){
   p[i] = 1/(alpha+(n-1));
  }

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
    new_utheta[0] = thetaD(y(ind,spanall),lambda_struct);
    thetaStar_vector[j] = new_utheta[0];
  }
  
  probs[nunique] = alpha/(alpha+n+0.0);
  indp = rmultinomF(probs);
  probssize = probs.size();
  if(indp == probssize) {
    outGD = GD(lambda_struct);
    thetaNp10_struct.mu = outGD.mu;
    thetaNp10_struct.rooti = outGD.rooti;
    thetaNp1_vector[0] = thetaNp10_struct;
  } else {
    outGD = thetaStar_vector[indp-1];
    thetaNp10_struct.mu = outGD.mu;
    thetaNp10_struct.rooti = outGD.rooti;
    thetaNp1_vector[0] = thetaNp10_struct;
  }
    
  //draw alpha
  alpha = alphaD(priorAlpha_struct,nunique,gridsize);

  //draw lambda
  lambda_struct = lambdaD(lambda_struct,thetaStar_vector,lambda_hyper["alim"],lambda_hyper["nulim"],lambda_hyper["vlim"],gridsize);

  DPOut out_struct;
    out_struct.indic = indic;
    out_struct.thetaStar_vector = thetaStar_vector;
    out_struct.thetaNp1_vector = thetaNp1_vector;
    out_struct.alpha = alpha;
    out_struct.Istar = nunique;
    out_struct.lambda_struct = lambda_struct;
    
  return(out_struct);
}

//RCPP SECTION----
//[[Rcpp::export]]
List rivDP_rcpp_loop(int R, int keep, int nprint,
                     int dimd, vec const& mbg, mat const& Abg, vec const& md, mat const& Ad,
                     vec const& y, bool isgamma, mat const& z, vec const& x, mat const& w, vec delta,
                     List const& PrioralphaList, int gridsize, bool SCALE, int maxuniq, double scalex, double scaley,
                     List const& lambda_hyper,double BayesmConstantA, int BayesmConstantnu){

  // Wayne Taylor 3/14/2015

  int n = y.size();
  int dimg = 1;
  if(isgamma) dimg = w.n_cols;

  //variable type initializaion
  int Istar, bgsize, mkeep;
  double beta;
  vec gammaVec, q0v, bg;
  mat errMat, wEmpty, V;
  wEmpty.reset(); //enforce 0 elements
  ytxtxtd out_struct;
  
  //initialize indicator vector, thetaStar, ncomp, alpha
  ivec indic = ones<ivec>(n);

  std::vector<murooti> thetaStar_vector(1), thetaNp1_vector(1);
  murooti thetaNp10_struct, thetaStar0_struct;
    thetaStar0_struct.mu = zeros<vec>(2);
    thetaStar0_struct.rooti = eye(2,2);
  thetaStar_vector[0] = thetaStar0_struct;
  
  //Initialize lambda
  lambda lambda_struct;
    lambda_struct.mubar = zeros<vec>(2);
    lambda_struct.Amu = BayesmConstantA;
    lambda_struct.nu = BayesmConstantnu;
    lambda_struct.V = lambda_struct.nu*eye(2,2);  
    
  //convert Prioralpha from List to struct
  priorAlpha priorAlpha_struct;
    priorAlpha_struct.power = PrioralphaList["power"];
    priorAlpha_struct.alphamin = PrioralphaList["alphamin"];
    priorAlpha_struct.alphamax = PrioralphaList["alphamax"];
    priorAlpha_struct.n = PrioralphaList["n"];  
  
  int ncomp = 1;
  
  double alpha = 1.0;
  
  //allocate space for draws
  mat deltadraw = zeros<mat>(R/keep,dimd);
  vec betadraw = zeros<vec>(R/keep);
  vec alphadraw = zeros<vec>(R/keep);
  vec Istardraw = zeros<vec>(R/keep);
  mat gammadraw = zeros<mat>(R/keep,dimg);
  List thetaNp1draw(R/keep);
  vec nudraw = zeros<vec>(R/keep);
  vec vdraw = zeros<vec>(R/keep);
  vec adraw = zeros<vec>(R/keep);

  if(nprint>0) startMcmcTimer();

  for(int rep = 0; rep < R; rep++) {
    
    //draw beta and gamma
    if(isgamma){
      out_struct = get_ytxt(y,z,delta,x,w,ncomp,indic,thetaStar_vector);
    } else {
      out_struct = get_ytxt(y,z,delta,x,wEmpty,ncomp,indic,thetaStar_vector);
    }
    
    bg = breg(out_struct.yt,out_struct.xt,mbg,Abg); 
    beta = bg[0];
    bgsize = bg.size()-1;
    if(isgamma) gammaVec = bg.subvec(1,bgsize);

    //draw delta
    if(isgamma){
      out_struct=get_ytxtd(y,z,beta,gammaVec,x,w,ncomp,indic,thetaStar_vector,dimd);
    } else {
      out_struct=get_ytxtd(y,z,beta,gammaVec,x,wEmpty,ncomp,indic,thetaStar_vector,dimd);
    }
  
    delta = breg(out_struct.yt,out_struct.xtd,md,Ad);
        
    //DP process stuff- theta | lambda
    if(isgamma) {
      errMat = join_rows(x-z*delta,y-beta*x-w*gammaVec);
    } else {
      errMat = join_rows(x-z*delta,y-beta*x);
    }
    
    q0v = q0(errMat,lambda_struct);
    
    DPOut DPout_struct = rthetaDP(maxuniq,alpha,lambda_struct,priorAlpha_struct,thetaStar_vector,indic,q0v,errMat,gridsize,lambda_hyper);
    
    indic = DPout_struct.indic;
    thetaStar_vector = DPout_struct.thetaStar_vector;
    alpha = DPout_struct.alpha;
    Istar = DPout_struct.Istar;
    thetaNp1_vector = DPout_struct.thetaNp1_vector;
    thetaNp10_struct = thetaNp1_vector[0];
    ncomp=thetaStar_vector.size();
    lambda_struct = DPout_struct.lambda_struct;
    
    if (nprint>0) if((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      deltadraw(mkeep-1,span::all) = trans(delta);
      betadraw[mkeep-1] = beta;
      alphadraw[mkeep-1] = alpha;
      Istardraw[mkeep-1] = Istar;
      if(isgamma) gammadraw(mkeep-1,span::all) = trans(gammaVec);
      //We need to convert from to NumericVector so that the nmix plotting works properly (it does not work for an nx1 matrix)
      thetaNp1draw[mkeep-1] = List::create(List::create(Named("mu") = NumericVector(thetaNp10_struct.mu.begin(),thetaNp10_struct.mu.end()),Named("rooti") = thetaNp10_struct.rooti));
      adraw[mkeep-1] = lambda_struct.Amu;
      nudraw[mkeep-1] = lambda_struct.nu;
      V = lambda_struct.V;
      vdraw[mkeep-1] = V(0,0)/(lambda_struct.nu+0.0);
    }
  }
  
  //rescale
  if(SCALE){
    deltadraw=deltadraw*scalex;
    betadraw=betadraw*scaley/scalex;
    if(isgamma) gammadraw=gammadraw*scaley;
  }  
  
  if (nprint>0) endMcmcTimer();
  
  return List::create(
    Named("deltadraw") = deltadraw,
    Named("betadraw") = betadraw,
    Named("alphadraw") = alphadraw,
    Named("Istardraw") = Istardraw,
    Named("gammadraw") = gammadraw,
    Named("thetaNp1draw") = thetaNp1draw,
    Named("adraw") = adraw,
    Named("nudraw") = nudraw,
    Named("vdraw") = vdraw);
}
