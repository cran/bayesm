#include "bayesm.h"
 
//SUPPORT FUNCTIONS SPECIFIC TO MAIN FUNCTION--------------------------------------------------------------------------------------
mat r2Sigma(vec const& r, int K){
//
// Keunwoo Kim 10/28/2014
//
// Purpose: 
//      convert r (vector) into Sigma (matrix)
//
// Arguments:
//      r : K*(K+1)/2 length vector
//      K : number of parameters (=nrow(Sigma))
//
// Output: 
//      Sigma (K by K matrix)
//
  int k, i, j;
  mat L = zeros<mat>(K, K);
  L.diag() = exp(r(span(0,K-1)));
  k = 0;
  for (i=0; i<K-1; i++){
	  for (j=i+1; j<K; j++){
		  L(j,i) = r[K+k];
		  k = k + 1;
	  }
  }
  return (L*trans(L));
}

double logJacob(mat const& choiceProb, int J){
//
// Keunwoo Kim 10/28/2014
//
// Purpose: 
//      compute log(det(Jacobian)) of mapping from ubobserved shock to share
//      (change-of-variables)
//
// Arguments:
//      choiceProb: T*J by H
//      J: number of alternatives (without outside option)
//
// Output: 
//      log(det(Jacobian)) of mapping from ubobserved shock to share
//
  int t;
  mat blockMat;
  double detblockMat;

  int H = choiceProb.n_cols;
  int T = choiceProb.n_rows / J;

  mat onesJJ = ones<mat>(J, J);  
  // equivalent to struc = kron(eye(T, T), onesJJ)
  mat struc = zeros<mat>(J*T,J*T);
  for (t=0; t<T; t++){
    struc(span(t*J,t*J+J-1),span(t*J,t*J+J-1)) = onesJJ;
  } 
  
  struc = struc - eye<mat>(T*J, T*J);
  mat offDiag = -choiceProb*trans(choiceProb)/H;
  mat Jac = struc%offDiag; 
  Jac.diag() = sum(choiceProb%(1-choiceProb), 1)/H;

  double sumlogJacob = 0;
  for (t=0; t<T; t++){
    blockMat = Jac(span(t*J, (t+1)*J-1), span(t*J, (t+1)*J-1)); 
    detblockMat = det(blockMat);
    // abs cannot be used for a scalor.
    sumlogJacob = sumlogJacob + log(sqrt(detblockMat*detblockMat));      
  }
  
  return (-sumlogJacob);  
}

mat share2mu(mat const& Sigma, mat const& X, mat const& v, vec const& share, int J, double tol){
//
// Keunwoo Kim 10/28/2014
//
// Purpose: 
//      contraction mapping (BLP)
//
// Arguments:
//      Sigma: var-cov matrix of random coefficients
//      X: J*T by K      
//      share: observed share (length J*T)
//      v: random draws from standard normal (K by H)
//      J: number of alternatives (without outside option)
//      tol: convergence tolerance for the contraction mapping
//
// Output: 
//      a matrix of mean utility (first column) and 
//      individual choice probabilities (second~last column)
//
  int t;
  mat expU, temp1, expSum, choiceProb;
  vec share_hat;  
  
  int H = v.n_cols;
  int T = X.n_rows/J;
  mat temp2(T*J,H);
  
  // T*J by H
  mat u = X*(trans(chol(Sigma))*v);
  int iter = 0;
  vec mu0 = ones<vec>(J*T);
  vec mu1 = mu0/2;
  
  //relative increasement
  vec rel = (mu1 - mu0)/mu0;
  double max_rel = max(abs(rel));
  while (max_rel > tol){
	  mu0 = mu1;
	  expU = exp(u + mu0*ones<mat>(1,H));

	  temp1 = reshape(expU, J, T*H);
	  expSum = 1 + sum(temp1, 0);
	  expSum = reshape(expSum, T, H);
	  // equivalent to expSum = kron(expSum, ones<vec>(J));    
    for (t=0; t<T; t++){
      temp2(span(t*J, t*J+J-1), span::all) = ones<vec>(J)*expSum(t, span::all);
    }
    expSum = temp2;
	  choiceProb = expU/expSum;
	  share_hat = sum(choiceProb, 1)/H;

	  mu1 = mu0 + log(share/share_hat);
	  iter = iter + 1;
	  rel = (mu0 - mu1)/mu0;
	  max_rel = max(abs(rel));
  }
  mat rtn = zeros(J*T, H+1);
  rtn(span::all,0) = mu1;
  rtn(span::all,span(1,H)) = choiceProb;

  return (rtn);
}

List rivDraw(vec const& mu, vec const& Xend, mat const& z, mat const& Xexo, vec const& theta_hat, mat const& A, 
                  vec const& deltabar, mat const& Ad, mat const& V, int nu, vec const& delta_old, mat const& Omega_old){
//
// Keunwoo Kim 05/21/2015
//
// Purpose: draw from posterior for linear I.V. model
//
// Arguments:
//        mu is vector of obs on lhs var in structural equation
//        Xend is "endogenous" var in structural eqn
//        Xexo is matrix of obs on "exogenous" vars in the structural eqn
//        z is matrix of obs on instruments
//
//        deltabar is prior mean of delta
//        Ad is prior prec
//        theta_hat is prior mean vector for theta2,theta1
//        A is prior prec of same
//        nu,V parms for IW on Omega
//
//        delta_old is the starting value from the previous chain
//        Omega_old is the starting value from the previous chain
//
// Output: list of draws of delta,thetabar,Omega
// 
// Model:
//    Xend=z'delta + e1
//    mu=thetabar1*Xend + Xexo'thetabar2 + e2
//        e1,e2 ~ N(0,Omega)
//
// Prior:
//   delta ~ N(deltabar,Ad^-1)
//   thetabar = vec(theta2,theta1) ~ N(theta_hat,A^-1)
//   Omega ~ IW(nu,V)
// 
  vec e1, ee2, bg, u, theta2;
  mat xt, Res, S, B, L, Li, z2, zt1, zt2, ucholinv, VSinv, mut;
  double sig,theta1;
  List out;
  int i;  

  int n = mu.size();
  int dimd = z.n_cols;
  int dimg = Xexo.n_cols;
  vec thetabar(dimg+1);

  mat C = eye(2,2);

  // set initial values
  mat Omega = Omega_old;
  vec delta = delta_old;

  mat xtd(2*n, dimd);  
  vec zvec = vectorise(trans(z));  
     
  //
  // draw beta,gamma
  //
  e1 = Xend - z*delta;
  ee2 = (Omega(0,1)/Omega(0,0)) * e1;
  sig = sqrt(Omega(1,1)-((Omega(0,1)*Omega(0,1))/Omega(0,0)));
  mut = (mu-ee2)/sig;
  xt = join_rows(Xend,Xexo)/sig;
  bg = breg(mut,xt,theta_hat,A);
  theta1 = bg[0];
  theta2 = bg(span(1,bg.size()-1));
    
  //
  // draw delta
  //
  C(1,0) = theta1;
  B = C*Omega*trans(C);
  L = trans(chol(B));
  Li = solve(trimatl(L),eye(2,2));
  u = mu - Xexo*theta2;
  mut = vectorise(Li * trans(join_rows(Xend,u)));
  z2 = trans(join_rows(zvec, theta1*zvec));
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
  delta = breg(mut,xtd,deltabar,Ad);
    
  //
  // draw Sigma
  //
  Res = join_rows(Xend-z*delta, mu-theta1*Xend-Xexo*theta2);
  S = trans(Res)*Res;
    
  // compute the inverse of V+S
  ucholinv = solve(trimatu(chol(V+S)), eye(2,2));
  VSinv = ucholinv*trans(ucholinv);
    
  out = rwishart(nu+n, VSinv);
  Omega = as<mat>(out["IW"]);
  
  thetabar(span(0,dimg-1)) = theta2;
  thetabar[dimg] = theta1;
  
  return List::create(
      Named("deltadraw") = delta,
      Named("thetabardraw") = thetabar,      
      Named("Omegadraw") = Omega
  );   
}

//MAIN FUNCTION---------------------------------------------------------------------------------------
// [[Rcpp::export]]
List bayesBLP_rcpp_loop(bool IV, mat const& X, mat const& Z, vec const& share, 
                        int J, int T, mat const& v, int R,
                        vec const& sigmasqR, 
                        mat const& A, vec const& theta_hat, 
                        vec const& deltabar, mat const& Ad,
                        int nu0, double s0_sq, mat const& VOmega, 
                        double ssq, mat const& cand_cov, 
                        vec const& theta_bar_initial, vec const& r_initial, 
                        double tau_sq_initial, mat const& Omega_initial, vec const& delta_initial,
                        double tol, int keep, int nprint){
//
// Keunwoo Kim 05/21/2015
//
// Purpose: 
//      draw theta_bar and Sigma via hybrid Gibbs sampler (Jiang, Manchanda, and Rossi, 2009)
//
// Arguments:
//    Observation
//      IV: whether to use instrumental variable (TRUE or FALSE)
//      X: J*T by H (If IV is TRUE, the last column is endogenous variable.)
//      z: instrumental variables (If IV is FALSE, it is not used.)
//      share: vector of length J*T
//
//    Dimension
//      J: number of alternatives
//      T: number of time
//      R: number of Gibbs sampling
//
//    Prior
//      sigmasqR
//      theta_hat
//      A
//      deltabar (used when IV is TRUE)
//      Ad (used when IV is TRUE)
//      nu0
//      s0_sq (used when IV is FALSE)
//      VOmega (used when IV is TRUE)
//
//    Metropolis-Hastings
//      ssq: scaling parameter
//      cand_cov: var-cov matrix of random walk
//
//    Initial values
//      theta_bar_initial
//      r_initial
//      tau_sq_initial (used when IV is FALSE)
//      Omega_initial (used when IV is TRUE)
//      delta_initial (used when IV is TRUE)
//
//    Contraction mapping
//      tol: convergence tolerance for the contraction mapping
//      v: draws used for Monte-Carlo integration
//
// Output:
//      a List of theta_bar, r (Sigma), tau_sq, Omega, and delta  draws
//      number of acceptance and loglikelihood
//
// Model & Prior: 
//      shown in the below comments.

  int nu1, mkeep, I, jt;
  mat prob_t, Sigma_new, b, S, Sigma, Sigma_inv, rel, expU, share_hat, choiceProb, expSum, L, ucholinv, XXAinv, out_cont,
      Xexo, Xend, Omega_all, delta_all, zetaeta_old, zetaeta_new, rootiOmega;
  vec r_new, mu_new, theta_tilde, z, mu, err, mu0, mu1, eta_new, eta_old, tau_sq_all, zeta;
  double alpha, ll_new, ll_old, sumLogJaco_new, prior_new, prior_old, s1_sq, acceptrate;
  List ivout;
  
  double pi = M_PI;
  int K = theta_hat.size();  
  
  if (IV==TRUE){
    Xexo = X(span::all, span(0,K-2));
    Xend = X(span::all, K-1);
    I = Z.n_cols;
  }
  
  // number of MC integration draws
  int H = v.n_cols;

  // Allocate matrix for draws to be stored during MCMC
  if (IV==TRUE){
    Omega_all = zeros<mat>(4,R/keep);
    delta_all = zeros<mat>(I,R/keep);
  }else{
    tau_sq_all = zeros<vec>(R/keep);
  }
  mat theta_bar_all = zeros<mat>(K,R/keep);
  mat r_all = zeros<mat>(K*(K+1)/2,R/keep);
  mat Sigma_all = zeros<mat>(K*K,R/keep);  
  vec ll_all = zeros<vec>(R/keep);

  // list to be returned to R  
  List rtn;

  // initial values
  vec theta_bar = theta_bar_initial;
  mat Omega = Omega_initial;
  vec delta = delta_initial;
  vec r_old = r_initial;
  double tau_sq = tau_sq_initial;  
  mat Sigma_old = r2Sigma(r_old, K);

  //===================================================================
  // get initial mu and sumLogJaco: Contraction Mapping
  //===================================================================
  // convert shares into mu
  out_cont = share2mu(Sigma_old, X, v, share, J, tol);
  mu = out_cont(span::all,0);
  choiceProb = out_cont(span::all,span(1,H));

  // Jacobian
  double sumLogJaco_old = logJacob(choiceProb, J);
  vec mu_old = mu;

  //===================================================================
  // Start MCMC
  //===================================================================
  if (nprint>0) startMcmcTimer();
  double n_accept = 0.0;
  for (int rep=0; rep<R; rep++){
	  //========================================================================
	  // STEP 1
	  // Draw r (for Sigma): Metropolis Hasting
	  // r_new = r_old + N(0, ssq*cand_cov)
	  // Prior:
	  // r ~ N(0, diag(sigmasqR)), that is, independent prior
	  //========================================================================
	  // get candidate
	  r_new = r_old + trans(chol(ssq*cand_cov))*randn<vec>(K*(K+1)/2);
	  Sigma_new = r2Sigma(r_new, K);
	  // convert share into mu_new
    out_cont = share2mu(Sigma_new, X, v, share, J, tol);
    mu_new = out_cont(span::all,0);
    choiceProb = out_cont(span::all,span(1,H));
	  // get eta_new
	  eta_new = mu_new - X*theta_bar;
	  // get eta_old
	  eta_old = mu_old - X*theta_bar;
    
    if (IV==TRUE){
      // get zeta
      zeta = Xend - Z*delta;
      // get ll_old
      zetaeta_old = join_rows(zeta, eta_old);
      rootiOmega = solve(trimatu(chol(Omega)), eye(2,2));
      ll_old = 0;
      for (jt=0; jt<J*T; jt++){
        ll_old = ll_old + lndMvn(vectorise(zetaeta_old(jt, span::all)), 
                          zeros<vec>(2), 
                          rootiOmega);
      }    
      ll_old = ll_old + sumLogJaco_old;
      // get ll_new    
      zetaeta_new = join_rows(zeta, eta_new);
	    sumLogJaco_new = logJacob(choiceProb, J);
      ll_new = 0;
      for (jt=0; jt<J*T; jt++){
        ll_new = ll_new + lndMvn(vectorise(zetaeta_new(jt, span::all)), 
                          zeros<vec>(2), 
                          rootiOmega);
      }    
      ll_new = ll_new + sumLogJaco_new;
    }else{
      // get ll_old
	    ll_old = sum(log((1/sqrt(2*pi*tau_sq)) * exp(-(eta_old%eta_old)/(2*tau_sq)))) + sumLogJaco_old;
	    // get ll_new	  
	    sumLogJaco_new = logJacob(choiceProb, J);
	    ll_new = sum(log((1/sqrt(2*pi*tau_sq)) * exp(-(eta_new%eta_new)/(2*tau_sq)))) + sumLogJaco_new;
    }
    
	  // priors
	  prior_new = sum(log((1/sqrt(2*pi*sigmasqR)) % exp(-(r_new%r_new)/(2*sigmasqR))));
	  prior_old = sum(log((1/sqrt(2*pi*sigmasqR)) % exp(-(r_old%r_old)/(2*sigmasqR))));

	  alpha = exp(ll_new + prior_new - ll_old - prior_old);
	  if (alpha>1) {alpha = 1;}

	  if (runif(1)[0]<=alpha) {
		  r_old = r_new;
		  Sigma_old = Sigma_new;
		  mu_old = mu_new;
		  sumLogJaco_old = sumLogJaco_new;
		  n_accept = n_accept + 1; 
	  }	
	  //========================================================================	
	  // STEP 2
	  // Draw theta_bar & tau^2 (or Omega & delta): Gibbs Sampler
    // mu = X*theta_bar + eta, eta~N(0,tau_sq)
    // (For IV case, see the comments in rivDraw above.)
	  // Prior:
	  // 1. theta_bar ~ N(theta_hat, A^-1)
	  // 2. tau_sq ~ nu0*s0_sq/chisq(nu0)
	  // Posterior:
	  // 1. theta_bar | tau_sq ~ N(theta_tilde, (X^t X/tau_sq + A)^-1)
	  // theta_tilde = (X^t X/tau_sq + A)^-1 * (tau_sq^-1*X^t mu + A*theta_hat)
	  // 2. tau_sq | theta_bar ~ nu1*s1_sq/chisq(nu1)
	  // nu1 = nu0 + n (n=J*T)
	  // s1_sq = [nu0*s0_sq + (mu-X^t theta_bar)^t (mu-X^t theta_bar)]/[nu0 + n]
	  //========================================================================    
    if (IV==TRUE){
      ivout = rivDraw(mu_old, Xend, Z, Xexo, theta_hat, A, 
                          deltabar, Ad, VOmega, nu0, delta, Omega);    
      delta = as<vec>(ivout["deltadraw"]);
      theta_bar = as<vec>(ivout["thetabardraw"]);
      Omega = as<mat>(ivout["Omegadraw"]);
    }else{    
      // compute the inverse of (trans(X)*X)/tau_sq + A
      ucholinv = solve(trimatu(chol((trans(X)*X)/tau_sq + A)), eye(K,K));
      XXAinv = ucholinv*trans(ucholinv);
 
      theta_tilde = XXAinv * (trans(X)*mu_old/tau_sq + A*theta_hat);
	    theta_bar = theta_tilde + ucholinv*vec(rnorm(K));

	    nu1 = nu0 + J*T;
	    err = mu_old - X*theta_bar;
	    s1_sq = (nu0*s0_sq + sum(err%err))/nu1;
	    z = vec(rnorm(nu1));	  
	    tau_sq = nu1*s1_sq/sum(z%z);
    }

	  //
    // print time to completion and draw # every nprint'th draw
    //
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
	  //========================================================================
	  // STEP 3
	  // Store Draws
	  //========================================================================
	  if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;      
      if (IV==TRUE){
        Omega_all(span::all,mkeep-1) = vectorise(Omega);      
        delta_all(span::all,mkeep-1) = delta;
      }else{
        tau_sq_all[mkeep-1] = tau_sq;
      }
	    theta_bar_all(span::all,mkeep-1) = theta_bar; 
      r_all(span::all,mkeep-1) = r_old;  
	    Sigma_all(span::all,mkeep-1) = vectorise(r2Sigma(r_old, K));	  
	    ll_all[mkeep-1] = ll_old;
	  }
  }
  acceptrate = n_accept/R;
  rtn["tausqdraw"] = tau_sq_all;
  rtn["Omegadraw"] = Omega_all;  
  rtn["deltadraw"] = delta_all;
  rtn["thetabardraw"] = theta_bar_all;
  rtn["rdraw"] = r_all;
  rtn["Sigmadraw"] = Sigma_all;  
  rtn["ll"] = ll_all;
  rtn["acceptrate"] = acceptrate;

  if (nprint>0) endMcmcTimer();
  return (rtn);
}
