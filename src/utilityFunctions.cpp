#include "bayesm.h"
 
//Used in rmvpGibbs and rmnpGibbs---------------------------------------------------------------------------------
vec condmom(vec const& x, vec const& mu, mat const& sigmai, int p, int j){
  
// Wayne Taylor 9/24/2014

//function to compute moments of x[j] | x[-j]
//output is a vec: the first element is the conditional mean
//                 the second element is the conditional sd

  vec out(2);
  int jm1 = j-1;
  int ind = p*jm1;
  
  double csigsq = 1./sigmai(ind+jm1);
  double m = 0.0;
  
  for(int i = 0; i<p; i++) if (i!=jm1) m += - csigsq*sigmai(ind+i)*(x[i]-mu[i]);
  
  out[0] = mu[jm1]+m;
  out[1] = sqrt(csigsq);
  
  return (out);
}

double rtrun1(double mu, double sigma,double trunpt, int above) {

// Wayne Taylor 9/8/2014
  
//function to draw truncated normal
//above=1 means from above b=trunpt, a=-inf
//above=0 means from below a=trunpt, b= +inf   
//modified by rossi 6/05 to check arg to qnorm

	double FA,FB,rnd,result,arg;
	if (above) {
		FA = 0.0; FB = R::pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
	} else {
		FB = 1.0; FA = R::pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
	}
	
  rnd = runif(1)[0]; //runif returns a NumericVector, so using [0] allows for conversion to double
	arg = rnd*(FB-FA)+FA;
	if(arg > .999999999) arg = .999999999;
	if(arg < .0000000001) arg = .0000000001;
	result = mu + sigma*R::qnorm(arg,0.0,1.0,1,0);

	return (result);
}

//Used in rhierLinearModel, rhierLinearMixture and rhierMnlRWMixture------------------------------------------------------
mat drawDelta(mat const& x,mat const& y,vec const& z,List const& comps,vec const& deltabar,mat const& Ad){

// Wayne Taylor 10/01/2014

// delta = vec(D)
//  given z and comps (z[i] gives component indicator for the ith observation, 
//   comps is a list of mu and rooti)
// y is n x p
// x is n x k
// y = xD' + U , rows of U are indep with covs Sigma_i given by z and comps

  int p = y.n_cols;
  int k = x.n_cols;
  int ncomp  = comps.length();
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
      
      List compsi = comps[compi];
      rowvec mui = as<rowvec>(compsi[0]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      mat rootii = trimatu(as<mat>(compsi[1])); //trimatu interprets the matrix as upper triangular
      yi.each_row() -= mui; //subtracts mui from each row of yi
      mat sigi = rootii*trans(rootii);
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

unireg runiregG(vec const& y, mat const& X, mat const& XpX, vec const& Xpy, double sigmasq, mat const& A, 
              vec const& Abetabar, int nu, double ssq) {

// Keunwoo Kim 09/16/2014

// Purpose: 
//  perform one Gibbs iteration for Univ Regression Model
//  only does one iteration so can be used in rhierLinearModel

// Model:
//  y = Xbeta + e  e ~N(0,sigmasq)
//  y is n x 1
//  X is n x k
//  beta is k x 1 vector of coefficients

// Prior:  
//  beta ~ N(betabar,A^-1)
//  sigmasq ~ (nu*ssq)/chisq_nu

  unireg out_struct;
  
  int n = y.size();
  int k = XpX.n_cols;
  
  //first draw beta | sigmasq
  mat IR = solve(trimatu(chol(XpX/sigmasq+A)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  vec btilde = (IR*trans(IR)) * (Xpy/sigmasq + Abetabar);
  vec beta = btilde + IR*vec(rnorm(k));
  
  //now draw sigmasq | beta
  double s = sum(square(y-X*beta));
  sigmasq = (s + nu*ssq)/rchisq(1,nu+n)[0]; //rchisq returns a vectorized object, so using [0] allows for the conversion to double
  
  out_struct.beta = beta;
  out_struct.sigmasq = sigmasq;  

  return (out_struct);
}

//Used in rnegbinRW and rhierNegbinRw-------------------------------------------------------------------------------------
double llnegbin(vec const& y, vec const& lambda, double alpha, bool constant){

// Keunwoo Kim 11/02/2014

// Computes the log-likelihood

// Arguments
//      y - a vector of observation
//      lambda - a vector of mean parameter (=exp(X*beta))
//      alpha - dispersion parameter
//      constant - TRUE(FALSE) if it computes (un)normalized log-likeihood

// PMF
//      pmf(y) = (y+alpha-1)Choose(y) * p^alpha * (1-p)^y
//      (y+alpha-1)Choose(y) = (alpha)*(alpha+1)*...*(alpha+y-1) / y! when y>=1 (0 when y=0)

  int i;
  int nobs = y.size();  
  vec prob = alpha/(alpha+lambda);    
  vec logp(nobs);
  if (constant){
    // normalized log-likelihood
    for (i=0; i<nobs; i++){
      // the fourth argument "1" indicates log-density
      logp[i] = R::dnbinom(y[i], alpha, prob[i], 1);
    }    
  }else{
    // unnormalized log-likelihood
    logp = sum(alpha*log(prob) + y % log(1-prob)); //% does element-wise multiplication
  }
  return (sum(logp));
}

double lpostbeta(double alpha, vec const& beta, mat const& X, vec const& y, vec const& betabar, mat const& rootA){

// Keunwoo Kim 11/02/2014

// Computes log posterior for beta | alpha

// Arguments
//        alpha - dispersion parameter of negative-binomial
//        beta - parameter of our interests
//        X, y - observation from data
//        betabar - mean of beta prior
//        rootA - t(rootA)%*%rootA = A (A^-1 is var-cov matrix of beta prior)

// Prior
//        beta ~ N(betabar, A^-1)

  vec lambda = exp(X*beta);
  double ll = llnegbin(y, lambda, alpha, FALSE);

  // unormalized prior
  vec z = rootA*(beta-betabar);
  double lprior = - 0.5*sum(z%z);
  
  return (ll+lprior);
}

double lpostalpha(double alpha, vec const& beta, mat const& X, vec const& y, double a, double b){

// Keunwoo Kim 11/02/2014

// Computes log posterior for alpha | beta

// Arguments
//        alpha - dispersion parameter of negative-binomial
//        beta - parameter of our interests
//        X, y - observation from data
//        a,b - parameters for Gamma distribution, alpha prior

// Prior
//        alpha ~ Gamma(a,b)
//        pdf(alpha) = b^a / Gamma(a) * alpha^(a-1) * e^(b*alpha)

  vec lambda = exp(X*beta);
  double ll = llnegbin(y, lambda, alpha, TRUE);
  // unormalized prior
  double lprior = (a-1)*log(alpha) - b*alpha;  
  
  return (ll+lprior);
}

//Used in rbprobitGibbs and rordprobitGibbs-----------------------------------------------------------------------
vec breg1(mat const& root, mat const& X, vec const& y, vec const& Abetabar) {

// Keunwoo Kim 06/20/2014

// Purpose: draw from posterior for linear regression, sigmasq=1.0

// Arguments:
//  root = chol((X'X+A)^-1)
//  Abetabar = A*betabar

// Output: draw from posterior

// Model: y = Xbeta + e  e ~ N(0,I)

// Prior: beta ~ N(betabar,A^-1)

  mat cov = trans(root)*root;  
    
  return (cov*(trans(X)*y+Abetabar) + trans(root)*vec(rnorm(root.n_cols)));
}

vec rtrunVec(vec const& mu,vec const& sigma, vec const& a, vec const& b){
  
// Keunwoo Kim 06/20/2014  

//function to draw from univariate truncated norm
//a is vector of lower bounds for truncation
//b is vector of upper bounds for truncation

  int n = mu.size();
  vec FA(n);
  vec FB(n);
  vec out(n);
  for (int i=0; i<n; i++) {
    FA[i] = R::pnorm((a[i]-mu[i])/sigma[i],0,1,1,0);
    FB[i] = R::pnorm((b[i]-mu[i])/sigma[i],0,1,1,0);
    out[i] = mu[i]+sigma[i]*R::qnorm(R::runif(0,1)*(FB[i]-FA[i])+FA[i],0,1,1,0);
  }

  return(out);
}

//Used in rhierMnlDP and rhierMnlRwMixture------------------------------------------------------------------------
mnlMetropOnceOut mnlMetropOnce(vec const& y, mat const& X, vec const& oldbeta, 
                                                 double oldll,double s, mat const& incroot, 
                                                 vec const& betabar, mat const& rootpi){ 
// Wayne Taylor 10/01/2014

// function to execute rw metropolis for the MNL
// y is n vector with element = 1,...,j indicating which alt chosen
// X is nj x k matrix of xvalues for each of j alt on each of n occasions
// RW increments are N(0,s^2*t(inc.root)%*%inc.root)
// prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
//  inc.root, rootpi are upper triangular
//  this means that we are using the UL decomp of Sigma^-1 for prior 
// oldbeta is the current


mnlMetropOnceOut metropout_struct;

double unif;
vec betadraw, alphaminv;

int stay = 0;
vec betac = oldbeta + s*trans(incroot)*as<vec>(rnorm(X.n_cols));
double cll = llmnl(betac,y,X);
double clpost = cll+lndMvn(betac,betabar,rootpi);
double ldiff = clpost-oldll-lndMvn(oldbeta,betabar,rootpi);
alphaminv << 1 << exp(ldiff);
double alpha = min(alphaminv);

     if(alpha < 1) {
       unif = runif(1)[0]; //runif returns a NumericVector, so using [0] allows for conversion to double
      } else { 
        unif=0;}
     if (unif <= alpha) {
       betadraw = betac;
       oldll = cll;
      } else {
        betadraw = oldbeta;
        stay = 1;
      }

metropout_struct.betadraw = betadraw;
metropout_struct.stay = stay;  
metropout_struct.oldll = oldll;

return (metropout_struct);
}

//Used in rDPGibbs, rhierMnlDP, rivDP-----------------------------------------------------------------------------
int rmultinomF(vec const& p){
  
// Wayne Taylor 1/28/2015

  vec csp = cumsum(p);
  double rnd = runif(1)[0]; //runif returns a NumericVector, so using [0] allows for conversion to double
  int res = 0;
  int psize = p.size();
  
  for(int i = 0; i < psize; i++){
    if(rnd > csp[i]) res = res+1;
  }
  
  return(res+1);
}

mat yden(std::vector<murooti> const& thetaStar_vector, mat const& y){

// Wayne Taylor 2/4/2015
  
// function to compute f(y | theta) 
// computes f for all values of theta in theta list of lists
      
// arguments:
//  thetaStar is a list of lists.  thetaStar[[i]] is a list with components, mu, rooti
//  y |theta[[i]] ~ N(mu,(rooti %*% t(rooti))^-1)  rooti is inverse of Chol root of Sigma

// output:
//  length(thetaStar) x n array of values of f(y[j,]|thetaStar[[i]]
  
  int nunique = thetaStar_vector.size();
  int n = y.n_rows;
  int k = y.n_cols;
  mat ydenmat = zeros<mat>(nunique,n);
  
  vec mu;
  mat rooti, transy, quads;
  
  for(int i = 0; i < nunique; i++){
    //now compute vectorized version of lndMvn 
    //compute y_i'RIRI'y_i for all i
        
    mu = thetaStar_vector[i].mu;
    rooti = thetaStar_vector[i].rooti;
  
    transy = trans(y);
    transy.each_col() -= mu; //column-wise subtraction
    
    quads = sum(square(trans(rooti) * transy),0); //same as colSums
    ydenmat(i,span::all) = exp(-(k/2.0)*log(2*M_PI) + sum(log(rooti.diag())) - .5*quads);
  }
  
  return(ydenmat);
}

ivec numcomp(ivec const& indic, int k){

// Wayne Taylor 1/28/2015
  
  //find the number of times each of k integers is in the vector indic
  ivec ncomp(k);
  
  for(int comp = 0; comp < k; comp++){
    ncomp[comp]=sum(indic == (comp+1));
  }
  
  return(ncomp);
}

murooti thetaD(mat const& y, lambda const& lambda_struct){

// Wayne Taylor 2/4/2015
  
// function to draw from posterior of theta given data y and base prior G0(lambda)
      
// here y ~ N(mu,Sigma)
// theta = list(mu=mu,rooti=chol(Sigma)^-1)
// mu|Sigma ~ N(mubar,Sigma (x) Amu-1)
// Sigma ~ IW(nu,V)
      
// arguments: 
//  y is n x k matrix of obs
//  lambda is list(mubar,Amu,nu,V)

// output:
//  one draw of theta, list(mu,rooti)
//  Sigma=inv(rooti)%*%t(inv(rooti))
      
// note: we assume that y is a matrix. if there is only one obs, y is a 1 x k matrix

  mat X = ones<mat>(y.n_rows,1);
  mat A(1,1); A.fill(lambda_struct.Amu);
  
  List rout = rmultireg(y,X,trans(lambda_struct.mubar),A,lambda_struct.nu,lambda_struct.V);
  
  murooti out_struct;
    out_struct.mu = as<vec>(rout["B"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    out_struct.rooti = solve(chol(trimatu(as<mat>(rout["Sigma"]))),eye(y.n_cols,y.n_cols)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  return(out_struct);
}

thetaStarIndex thetaStarDraw(ivec indic, std::vector<murooti> thetaStar_vector, mat const& y, mat ydenmat, vec const& q0v, double alpha, 
                        lambda const& lambda_struct, int maxuniq) {
                          
// Wayne Taylor 2/4/2015
                               
// indic is n x 1 vector of indicator of which of thetaStar is assigned to each observation
// thetaStar is list of the current components (some of which may never be used)
// y is n x d matrix of observations
// ydenmat is maxuniq x n matrix to store density evaluations - we assume first 
// length(Thetastar) rows are filled in with density evals
// q0v is vector of bayes factors for new component and each observation
// alpha is DP process tightness prior
// lambda is list of priors for the base DP process measure
// maxuniq maximum number of mixture components
// yden is function to fill out an array
// thetaD is function to draw theta from posterior of theta given y and G0
   
  int n = indic.size();
  ivec ncomp, indicC;
  int k, inc, cntNonzero;
  std::vector<murooti> listofone_vector(1);
  std::vector<murooti> thetaStarC_vector;
  
  //draw theta_i given theta_-i
  for(int i = 0; i<n; i++){
   k = thetaStar_vector.size();
   vec probs(k+1);
   probs[k] = q0v[i]*(alpha/(alpha+(n-1)));
   
   //same as to indicmi = indic[-i]
   ivec indicmi = zeros<ivec>(n-1);
   inc = 0;
   for(int j = 0; j<(n-1); j++){
     if(j == i) {inc = inc + 1;}
     indicmi[j] = indic[inc];
     inc = inc+1;
   }
   
   ncomp = numcomp(indicmi,k);
   
   for(int comp = 0; comp<k; comp++){
     probs[comp] = ydenmat(comp,i)*ncomp[comp]/(alpha+(n-1));
   }
   
   probs = probs/sum(probs);
   indic[i] = rmultinomF(probs);
  
   if(indic[i] == (k+1)){
     if((k+1) > maxuniq) {
        stop("max number of comps exceeded");
     } else {
      listofone_vector[0] = thetaD(y(i,span::all),lambda_struct);
      thetaStar_vector.push_back(listofone_vector[0]);
      ydenmat(k,span::all) = yden(listofone_vector,y);
    }}
  }
  
  //clean out thetaStar of any components which have zero observations associated with them
  //and re-write indic vector 
  k = thetaStar_vector.size();
  indicC = zeros<ivec>(n);
  ncomp = numcomp(indic,k);
  
  cntNonzero = 0;
  for(int comp = 0; comp<k; comp++){
   if(ncomp[comp] != 0){
     thetaStarC_vector.push_back(thetaStar_vector[comp]);
     cntNonzero=cntNonzero+1;
   for(int i = 0; i<n; i++){if(indic[i] == (comp+1)) indicC[i] = cntNonzero;} //same as indicC(indic==comp) = cntNonzero;
   }
  }

  thetaStarIndex out_struct;
    out_struct.indic = indicC;
    out_struct.thetaStar_vector = thetaStarC_vector;

  return(out_struct);
}

vec q0(mat const& y, lambda const& lambda_struct){
  
// Wayne Taylor 2/4/2015

// function to compute a vector of int f(y[i]|theta) p(theta|lambda)dlambda
// here p(theta|lambda) is G0 the base prior

// implemented for a multivariate normal data density and standard conjugate prior:
//  theta=list(mu,Sigma)
//  f(y|theta,eta) is N(mu,Sigma)
//  lambda=list(mubar,Amu,nu,V)
//    mu|Sigma ~ N(mubar,Sigma (x) Amu^-1)
//    Sigma ~ IW(nu,V)

// arguments:
//  Y is n x k matrix of observations
//  lambda=list(mubar,Amu,nu,V)
 
// output:
//  vector of q0 values for each obs (row of Y)

// p. rossi 12/05
//  here y is matrix of observations (each row is an obs)
  
  int k = y.n_cols;
  mat R = chol(lambda_struct.V);
  double logdetR = sum(log(R.diag()));
  double lnk1k2, constant;
  mat transy, m, vivi, lnq0v;
  
  if (k > 1) {
    vec km1(k-1); for(int i = 0; i < (k-1); i++) km1[i] = i+1; //vector of 1:k SEE SEQ_ALONG
    lnk1k2 = (k/2.0)*log(2.0)+log((lambda_struct.nu-k)/2)+lgamma((lambda_struct.nu-k)/2)-lgamma(lambda_struct.nu/2)+sum(log(lambda_struct.nu/2-km1/2));
  } else {
    lnk1k2 = (k/2.0)*log(2.0)+log((lambda_struct.nu-k)/2)+lgamma((lambda_struct.nu-k)/2)-lgamma(lambda_struct.nu/2);
  }
  
  constant = -(k/2.0)*log(2*M_PI)+(k/2.0)*log(lambda_struct.Amu/(1+lambda_struct.Amu)) + lnk1k2 + lambda_struct.nu*logdetR;

// note: here we are using the fact that |V + S_i | = |R|^2 (1 + v_i'v_i)
//  where v_i = sqrt(Amu/(1+Amu))*t(R^-1)*(y_i-mubar), R is chol(V)
//  and S_i = Amu/(1+Amu) * (y_i-mubar)(y_i-mubar)'
      
  transy = trans(y);
  transy.each_col() -= lambda_struct.mubar;
  
  m = sqrt(lambda_struct.Amu/(1+lambda_struct.Amu))*trans(solve(trimatu(R),eye(y.n_cols,y.n_cols)))*transy; //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  vivi = sum(square(m),0);
  
  lnq0v = constant - ((lambda_struct.nu+1)/2)*(2*logdetR+log(1+vivi));
  
  return(trans(exp(lnq0v)));
}

vec seq_rcpp(double from, double to, int len){

// Wayne Taylor 1/28/2015

// Same as R::seq()

  vec res(len);
  res[len-1] = to; res[0] = from; //note the order of these two statements is important, when gridsize = 1 res[0] will be rewritten to the correct number
  double increment = (res[len-1]-res[0])/(len-1);
  for(int i = 1; i<(len-1); i++) res[i] = res[i-1] + increment;
  return(res);
}

double alphaD(priorAlpha const& priorAlpha_struct, int Istar, int gridsize){

// Wayne Taylor 2/4/2015
  
// function to draw alpha using prior, p(alpha)= (1-(alpha-alphamin)/(alphamax-alphamin))**power
      
  //same as seq
  vec alpha = seq_rcpp(priorAlpha_struct.alphamin,priorAlpha_struct.alphamax-.000001,gridsize);
  
  vec lnprob(gridsize);
  for(int i = 0; i<gridsize; i++){
    lnprob[i] = Istar*log(alpha[i]) + lgamma(alpha[i]) - lgamma(priorAlpha_struct.n+alpha[i]) + priorAlpha_struct.power*log(1-(alpha[i]-priorAlpha_struct.alphamin)/(priorAlpha_struct.alphamax-priorAlpha_struct.alphamin));
  }
  
  lnprob = lnprob - median(lnprob);
  vec probs=exp(lnprob);
  probs=probs/sum(probs);
  
  return(alpha(rmultinomF(probs)-1));
}

murooti GD(lambda const& lambda_struct){
  
// Wayne Taylor 2/4/2015
      
// function to draw from prior for Multivariate Normal Model
      
// mu|Sigma ~ N(mubar,Sigma x Amu^-1)
// Sigma ~ IW(nu,V)

// note: we must insure that mu is a vector to use most efficient lndMvn routine

  int k = lambda_struct.mubar.size();
  
  List Rout = rwishart(lambda_struct.nu,solve(trimatu(lambda_struct.V),eye(k,k))); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  mat Sigma = as<mat>(Rout["IW"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
  mat root = chol(Sigma);
  mat draws = rnorm(k);
  mat mu = lambda_struct.mubar + (1/sqrt(lambda_struct.Amu))*trans(root)*draws;
  
  murooti out_struct;
    out_struct.mu = mu;
    out_struct.rooti = solve(trimatu(root),eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient

  return(out_struct);
}


lambda lambdaD(lambda const& lambda_struct, std::vector<murooti> const& thetaStar_vector, vec const& alim, vec const& nulim, vec const& vlim, int gridsize){

// Wayne Taylor 2/4/2015

// revision history
//  p. rossi 7/06
//  vectorized 1/07
//  changed 2/08 to paramaterize V matrix of IW prior to nu*v*I; then mode of Sigma=nu/(nu+2)vI
//    this means that we have a reparameterization to v* = nu*v

// function to draw (nu, v, a) using uniform priors

// theta_j=(mu_j,Sigma_j)  mu_j~N(0,Sigma_j/a)  Sigma_j~IW(nu,vI)
//  recall E[Sigma]= vI/(nu-dim-1)

  vec lnprob, probs, rowSumslgammaarg;
  int ind; //placeholder for matrix indexing
  murooti thetaStari_struct; mat rootii; vec mui;
  mat mout, rimu, arg, lgammaarg;
  double sumdiagriri, sumlogdiag, sumquads, adraw, nudraw, vdraw;

  murooti thetaStar0_struct = thetaStar_vector[0];
  int d = thetaStar0_struct.mu.size();
  int Istar = thetaStar_vector.size();
  
  vec aseq = seq_rcpp(alim[0],alim[1],gridsize);
  vec nuseq = d-1+exp(seq_rcpp(nulim[0],nulim[1],gridsize)); //log uniform grid
  vec vseq = seq_rcpp(vlim[0],vlim[1],gridsize);

// "brute" force approach would simply loop over the 
//  "observations" (theta_j) and use log of the appropriate densities.  To vectorize, we
// notice that the "data" comes via various statistics:
//  1. sum of log(diag(rooti_j)
//  2. sum of tr(V%*%rooti_j%*%t(rooti_j)) where V=vI_d
//  3. quadratic form t(mu_j-0)%*%rooti%*%t(rooti)%*%(mu_j-0)
// thus, we will compute these first.
// for documentation purposes, we leave brute force code in comment fields

// extract needed info from thetastar list
  
  //mout has the rootis in form: [t(rooti_1), t(rooti_2), ...,t(rooti_Istar)]
  mout = zeros<mat>(d,Istar*d);
  ind = 0;
  for(int i = 0; i < Istar; i++){
    thetaStari_struct = thetaStar_vector[i];
    rootii = thetaStari_struct.rooti;
    ind = i*d;
    mout.submat(0, ind,d-1,ind+d-1) = trans(rootii);
  }
  sumdiagriri = sum(sum(square(mout),0)); //sum_i trace(rooti_i*trans(rooti_i))

// now get diagonals of rooti
  sumlogdiag = 0.0;
  for(int i = 0; i < Istar; i++){
    ind = i*d;
    for(int j = 0; j < d; j++){
      sumlogdiag = sumlogdiag+log(mout(j,ind+j));
    }
  }

  //columns of rimu contain trans(rooti_i)*mu_i
  rimu = zeros<mat>(d,Istar);
  for(int i = 0; i < Istar; i++){
    thetaStari_struct = thetaStar_vector[i];
    mui = thetaStari_struct.mu;
    rootii = thetaStari_struct.rooti;
    rimu(span::all,i) = trans(rootii) * mui;
  }
  sumquads = sum(sum(square(rimu),0));

// draw a  (conditionally indep of nu,v given theta_j)
  lnprob = zeros<vec>(aseq.size());
// for(i in seq(along=aseq)){
// for(j in seq(along=thetastar)){
// lnprob[i]=lnprob[i]+lndMvn(thetastar[[j]]$mu,c(rep(0,d)),thetastar[[j]]$rooti*sqrt(aseq[i]))}
  lnprob = Istar*(-(d/2.0)*log(2*M_PI))-.5*aseq*sumquads+Istar*d*log(sqrt(aseq))+sumlogdiag;
  lnprob = lnprob-max(lnprob) + 200;
  probs = exp(lnprob);
  probs = probs/sum(probs);
  adraw = aseq[rmultinomF(probs)-1];

// draw nu given v

  lnprob = zeros<vec>(nuseq.size());
// for(i in seq(along=nuseq)){
// for(j in seq(along=thetastar)){
// Sigma_j=crossprod(backsolve(thetastar[[j]]$rooti,diag(d)))
// lnprob[i]=lnprob[i]+lndIWishart(nuseq[i],V,Sigma_j)}

  //same as arg = (nuseq+1-arg)/2.0;
  arg = zeros<mat>(gridsize,d);
  for(int i = 0; i < d; i++) {
    vec indvec(gridsize);
    indvec.fill(-(i+1)+1);
    arg(span::all,i) = indvec;
  }
  arg.each_col() += nuseq;
  arg = arg/2.0;

  lgammaarg = zeros<mat>(gridsize,d);
  for(int i = 0; i < gridsize; i++){
    for(int j = 0; j < d; j++){
      lgammaarg(i,j) = lgamma(arg(i,j));
  }}
  rowSumslgammaarg = sum(lgammaarg,1);
  
  lnprob = zeros<vec>(gridsize);
  for(int i = 0; i<gridsize; i++){
    lnprob[i] = -Istar*log(2.0)*d/2.0*nuseq[i] - Istar*rowSumslgammaarg[i] + Istar*d*log(sqrt(lambda_struct.V(0,0)))*nuseq[i] + sumlogdiag*nuseq[i];
  }
  
  lnprob = lnprob-max(lnprob)+200;
  probs = exp(lnprob);
  probs = probs/sum(probs);
  nudraw = nuseq[rmultinomF(probs)-1];

// draw v given nu 
      
  lnprob = zeros<vec>(vseq.size());
// for(i in seq(along=vseq)){
// V=vseq[i]*diag(d)
// for(j in seq(along=thetastar)){
// Sigma_j=crossprod(backsolve(thetastar[[j]]$rooti,diag(d)))
// lnprob[i]=lnprob[i]+lndIWishart(nudraw,V,Sigma_j)}
// lnprob=Istar*nudraw*d*log(sqrt(vseq))-.5*sumdiagriri*vseq
      
  lnprob = Istar*nudraw*d*log(sqrt(vseq*nudraw))-.5*sumdiagriri*vseq*nudraw;
  lnprob = lnprob-max(lnprob)+200;
  probs = exp(lnprob);
  probs = probs/sum(probs);
  vdraw = vseq[rmultinomF(probs)-1];

// put back into lambda
  lambda out_struct;
    out_struct.mubar = zeros<vec>(d);
    out_struct.Amu = adraw;
    out_struct.nu = nudraw;
    out_struct.V = nudraw*vdraw*eye(d,d);
  
  return(out_struct);
}

//Used in llnhlogit and simnhlogit---------------------------------------------------------------------------------
double root(double c1, double c2, double tol, int iterlim){

//function to find root of c1 - c2u = lnu
   
   int iter = 0;
   double uold = .1;
   double unew = .00001;
   
   while (iter <= iterlim && fabs(uold-unew) > tol){
     uold = unew;
     unew=uold + (uold*(c1 -c2*uold -  log(uold)))/(1. + c2*uold);
     if(unew < 1.0e-50) unew=1.0e-50;
     iter=iter+1;
   }
   
   return(unew);
}

//[[Rcpp::export]]
vec callroot(vec const& c1, vec const& c2, double tol, int iterlim){
  
  int n = c1.size();
  vec u = zeros<vec>(n);
  
  for(int i = 0; i<n; i++){
    u[i] = root(c1[i],c2[i],tol,iterlim);
  }
  
  return(u);
}
