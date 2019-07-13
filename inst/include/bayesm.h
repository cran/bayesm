#ifndef __BAYESM_H__
#define __BAYESM_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>

using namespace arma;
using namespace Rcpp;

//CUSTOM STRUCTS--------------------------------------------------------------------------------------------------
//Used in rhierLinearMixture, rhierLinearModel, rhierMnlDP, rhierMnlRwMixture, rhierNegbinRw, and rsurGibbs
struct moments{
  vec y;
  mat X;
  mat XpX;
  vec Xpy;
  mat hess;
};

//Used in rhierLinearMixture, rhierLinearModel, rhierMnlRWMixture, and utilityFunctions.cpp
struct unireg{
    vec beta;
    double sigmasq;
  };

//Used in rhierMnlDP, rhierMnlRwMixture, and utilityFunctions.cpp
struct mnlMetropOnceOut{
  vec betadraw;
  int stay;
  double oldll;
};  
  
//Used in rDPGibbs, rhierMnlDP, rivDP, and utilityFunctions.cpp
struct lambda{
    vec mubar;
    double Amu;
    double nu;
    mat V;
};

//Used in rDPGibbs, rhierMnlDP, rivDP, and utilityFunctions.cpp
struct priorAlpha{
  double power;
  double alphamin;
  double alphamax;
  int n;
};

//Used  in rDPGibbs, rhierMnlDP, rivDP, and utilityFunctions.cpp
struct murooti{
  vec mu;
  mat rooti;
};

//Used in rDPGibbs, rhierMnlDP, rivDP, and utilityFunctions.cpp
struct thetaStarIndex{
  ivec indic;
  std::vector<murooti> thetaStar_vector;
};

//Used in rhierMnlDP, rivDP
struct DPOut{
  ivec indic;
  std::vector<murooti> thetaStar_vector;
  std::vector<murooti> thetaNp1_vector;
  double alpha;
  int Istar;
  lambda lambda_struct;
};

//EXPOSED FUNCTIONS-----------------------------------------------------------------------------------------------
List rwishart(double nu, mat const& V);

List rmultireg(mat const& Y, mat const& X, mat const& Bbar, mat const& A, double nu, mat const& V);

vec rdirichlet(vec const& alpha);

double llmnl(vec const& beta, vec const& y, mat const& X);

mat lndIChisq(double nu, double ssq, mat const& X);

double lndMvst(vec const& x, double nu, vec const& mu, mat const& rooti, bool NORMC);

double lndMvn(vec const& x, vec const& mu, mat const& rooti);

double lndIWishart(double nu, mat const& V, mat const& IW);

vec rmvst(double nu, vec const& mu, mat const& root);

vec breg(vec const& y, mat const& X, vec const& betabar, mat const& A);

vec cgetC(double e, int k);

List rmixGibbs( mat const& y,  mat const& Bbar, mat const& A, double nu, mat const& V,  vec const& a, vec const& p,  vec const& z);
  //rmixGibbs contains the following support functions, which are called ONLY THROUGH rmixGibbs: drawCompsFromLabels, drawLabelsFromComps, and drawPFromLabels

//SUPPORT FUNCTIONS (contained in utilityFunctions.cpp and trunNorm.cpp)-----------------------------------------------------------
//Used in rmvpGibbs and rmnpGibbs
vec condmom(vec const& x, vec const& mu, mat const& sigmai, int p, int j);

//double rtrun1(double mu, double sigma,double trunpt, int above); <--NO LONGER USED

double trunNorm(double mu,double sig, double trunpt, int above);

//Used in rhierLinearModel, rhierLinearMixture and rhierMnlRWMixture
mat drawDelta(mat const& x,mat const& y,vec const& z,List const& comps,vec const& deltabar,mat const& Ad);

unireg runiregG(vec const& y, mat const& X, mat const& XpX, vec const& Xpy, double sigmasq, mat const& A, vec const& Abetabar, double nu, double ssq);

//Used in rnegbinRW and rhierNegbinRw
double llnegbin(vec const& y, vec const& lambda, double alpha, bool constant);

double lpostbeta(double alpha, vec const& beta, mat const& X, vec const& y, vec const& betabar, mat const& rootA);

double lpostalpha(double alpha, vec const& beta, mat const& X, vec const& y, double a, double b);

//Used in rbprobitGibbs (uses breg1 and trunNorm_vec) and rordprobitGibbs (uses breg1 and rtrunVec)
vec breg1(mat const& root, mat const& X, vec const& y, vec const& Abetabar);

vec rtrunVec(vec const& mu,vec const& sigma, vec const& a, vec const& b);

vec trunNorm_vec(vec const& mu, vec const& sig, vec const& trunpt, vec const& above);

//Used in rhierMnlDP and rhierMnlRwMixture
mnlMetropOnceOut mnlMetropOnce(vec const& y, mat const& X, vec const& oldbeta, double oldll,double s, mat const& incroot, vec const& betabar, mat const& rootpi);

//Used in rDPGibbs, rhierMnlDP, rivDP
int rmultinomF(vec const& p);

mat yden(std::vector<murooti> const& thetaStar, mat const& y);

ivec numcomp(ivec const& indic, int k);

murooti thetaD(mat const& y, lambda const& lambda_struct);

thetaStarIndex thetaStarDraw(ivec indic, std::vector<murooti> thetaStar_vector, mat const& y, mat ydenmat, vec const& q0v, double alpha, lambda const& lambda_struct, int maxuniq);

vec q0(mat const& y, lambda const& lambda_struct);

vec seq_rcpp(double from, double to, int len); //kept _rcpp due to conflict with base seq function

double alphaD(priorAlpha const& priorAlpha_struct, int Istar, int gridsize);

murooti GD(lambda const& lambda_struct);

lambda lambdaD(lambda const& lambda_struct, std::vector<murooti> const& thetaStar_vector, vec const& alim, vec const& nulim, vec const& vlim, int gridsize);

//FUNCTION TIMING (contained in functionTiming.cpp)---------------------------------------------------------------
void startMcmcTimer();
void infoMcmcTimer(int rep, int R);
void endMcmcTimer();

#endif
