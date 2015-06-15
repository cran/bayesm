#include "bayesm.h"
#include <RcppArmadilloExtensions/sample.h> //used for "sample" function
 
//SUPPORT FUNCTIONS SPECIFIC TO MAIN FUNCTION--------------------------------------------------------------------------------------
double ghk(mat const& L, vec const& a, vec const& b, int const& n, int const& dim){

//  Wayne Taylor 4/29/15

//routine to implement ghk with a region : a[i-1] <= x_i <= b[i-1]
//r mcculloch 8/04
//L is lower triangular root of Sigma random vector is assumed to have zero mean
//n is number of draws to use in GHK
//dim is the dimension of L
//modified 6/05 by rossi to check arg into qnorm
//converted to rcpp 5/15

  int i,j;
  
  NumericVector aa(1),bb(1),pa(1),pb(1),arg(1);
  double u,prod,mu;
  vec z(dim);
  
  double res=0.0;
  
  for(i=0;i<n;i++) {
      
    prod = 1.0;
      
    for(j=0;j<dim;j++) {
         
      mu = 0.0; 
      if(j>0) mu = as_scalar(L(j,span(0,j-1))*z(span(0,j-1))); //previously done via a loop for(k=0;k<j;k++) mu += L(k*dim+j)*z[k];
      
      aa[0] = (a[j]-mu)/L(j,j); //when using one element length NumericVectors, use [0] as much as possible
      bb[0] = (b[j]-mu)/L(j,j);
      
      pa[0] = pnorm(aa,0.0,1.0)[0];
      pb[0] = pnorm(bb,0.0,1.0)[0];
       
      prod *= pb[0]-pa[0];
      
      u = unif_rand(); //unif_rand() is slightly faster than runif() for single draws double
      
      arg[0] = u*pb[0]+(1.0-u)*pa[0];
      
      if(arg[0] > .999999999) arg[0]=.999999999;
      if(arg[0] < .0000000001) arg[0]=.0000000001;
      
      z[j] = qnorm(arg,0.0,1.0)[0];
    }
    
    res += prod;
  }
  
  res /= n; 
  
  return (res);
}

mat dy(mat y, mat const& x, vec const& c, vec const& mu, mat const& beta, vec const& s, vec const& tau, vec const& sigma){

// Wayne Taylor 4/29/15

  //Variable declaration
  double sigman,taun;
  rowvec yn;
  vec xn;
  
  int p = y.n_cols;
  int nobs = y.n_rows;
   
  //cm = conditional mean, cs = condtional standard deviation
  //u =uniform for truncated normal draw
  double cm,cs,u;
  
  // standardized truncation points (a,b)
  // cdf at truncation points (pa,pb)
  NumericVector a(1),b(1),pa(1),pb(1);
  double qout;

  //loop over coordinates of y - first by rows, then by columns
  for(int n = 0; n<nobs; n++){
    
    sigman = sigma[n];
    taun = tau[n];
    yn = trans(vectorise(y(n,span::all)));
    xn = vectorise(x(n,span::all));
    
    for(int i = 0; i<p; i++) {
      
      //compute conditonal mean and standard deviation
      cs = s[i]*sigman;
      cm = mu[i]+taun;
      
      for(int j=0;j<i;j++) cm += (beta(i*(p-1)+j))*(yn[j]-mu[j]-taun);
      for(int j=(i+1);j<p;j++) cm += (beta(i*(p-1)+j-1))*(yn[j]-mu[j]-taun);

      //draw truncated normal
      // y~N(cm,cs^2) I[c[x[i]-1],c[x[i])
      //a = (c[x(_,i)-1]-cm)/cs;  b = (c[x(_,i)]-cm)/cs;
      a[0] = ((c[xn[i]-1]-cm)/cs); //when using one element length NumericVectors, use [0] as much as possible
      b[0] = ((c[xn[i]]-cm)/cs);
      
      pa[0] = pnorm(a,0.0,1.0)[0];
      pb[0] = pnorm(b,0.0,1.0)[0];
      
      u = unif_rand(); //unif_rand() is slightly faster than runif() for single draws double
      
      qout = qnorm(u*pb + (1-u)*pa,0.0,1.0)[0];
      yn[i] = cm + cs*qout;
    }
    
    //put yn values back into y
    y(n,span::all) = yn;
  }
   
  return(y);
}

double rlpx(mat const& x, double e,int k, vec const& mu,vec const& tau,mat const& Sigma,vec const& sigma,int nd=500) {
  
//Wayne Taylor 4/29/15

  int n = x.n_rows;
  int p = x.n_cols;
  vec cc = cgetC(e,k);
  mat L = trans(chol(Sigma));
  vec lpv = zeros<vec>(n);
  double offset = p*log((double)k);

  vec a,b;
  double ghkres,lghkres;
  uvec xia(p),xib(p);
  mat Li;
  
  for(int i = 0; i<n; i++){
    Li = sigma[i]*L;
  
    for(int u = 0;u<p;u++){
      xia[u] = x(i,u)-1;
      xib[u] = x(i,u);
    }
    
    a = cc.elem(xia)-mu-tau[i];
    b = cc.elem(xib)-mu-tau[i];
    
    ghkres = ghk(Li,a,b,nd,L.n_rows);
    lghkres = trunc_log(ghkres); //natural log, truncated to avoid +/- infinity. Note on my machine it truncates to ~log(1e-308)
    lpv[i] = lghkres + offset;
  }
  
  return(sum(lpv));
}

List condd(mat const& Sigma) {
  
//Wayne Taylor 4/29/15
  
  int p = Sigma.n_rows;
  mat Si = solve(Sigma,eye(p,p));
  int cbetarows = p-1;
  mat cbeta = zeros<mat>(cbetarows,p);
  uvec ui(1),ind(p-1);
  int counter;

  uvec cbetaAllRow(cbetarows);
  for(int i = 0; i<cbetarows; i++) cbetaAllRow[i] = i;    

  for(int i = 0; i<p; i++){
    ui[0] = i;
    
    counter = 0;
    for(int j = 0; j<cbetarows; j++){
      if(j==i) counter = counter + 1;
      ind[j] = counter;
      counter = counter + 1;  
    }

    cbeta(cbetaAllRow,ui) = -Si(ind,ui)/as_scalar(Si(ui,ui));
  }
  
  return List::create(
    Named("beta") = cbeta,
    Named("s") = sqrt(1/Si.diag()));  
}

mat getS(mat const& Lam, int n, vec const& moms){
  
//Wayne Taylor 4/29/15
  
  mat S = zeros<mat>(2,2);
  
  S(0,0) = (n-1)*moms[2] + n*pow(moms[0],2);
  S(0,1) = (n-1)*moms[3] + n*moms[0]*(moms[1]-Lam(1,1));
  S(1,0) = S(0,1);
  S(1,1) = (n-1)*moms[4] + n*pow(moms[1]-Lam(1,1),2);

  return(S);
}

double llL(mat const& Lam, int n, mat const& S, mat const& V,int nu){
  
//Wayne Taylor 4/29/15  

  int d = Lam.n_cols;
  double dlam = Lam(0,0)*Lam(1,1)-pow(Lam(0,1),2);
  mat M = (S+V) *  solve(Lam,eye(d,d));
  double ll = -.5*(n+nu+3)*log(dlam) -.5*sum(M.diag());
  
  return(ll);
}

//MAIN FUNCTION------------------------------------------------------------------------------------
//[[Rcpp::export]]
List rscaleUsage_rcpp_loop(int k, mat const& x, int p, int n,
                           int R, int keep, int ndghk, int nprint,
                           mat y, vec mu, mat Sigma, vec tau, vec sigma, mat Lambda, double e,
                           bool domu, bool doSigma, bool dosigma, bool dotau, bool doLambda, bool doe,
                           int nu, mat const& V, mat const& mubar, mat const& Am,
                           vec const& gsigma, vec const& gl11,vec const& gl22, vec const& gl12,
                           int nuL, mat const& VL, vec const& ge){

// R.McCulloch, 12/04  code for scale usage R function (rScaleUsage)
//  changed to R error function, P. Rossi 05/12
//  converted to rcpp W. Taylor 04/15

  //variable declaration
  int mkeep, ng, ei, pi;
  double eprop, eold;
  double Ai, A, xtx, beta, s2, m, a, b, s, qr, llold, llprop, lrat, paccept;
  vec cc, xty, ete, pv, h, moms(5),rgl11, rgl12a, rgl12, rgl22, absege, minvec;
  uvec eiu;
  mat Res, S, yd, Si, Vmi, Rm, Ri, Vm, mm, onev, xx, ytemp, yy, eps, dat, temp, SS;
  List bs, rwout;
  
  rowvec onesp = ones<rowvec>(p);
  int nk = R/keep;
  int ndpost = nk*keep;
  
  mat drSigma = zeros<mat>(nk,pow(p,2.0));
  mat drmu = zeros<mat>(nk,p);
  mat drtau = zeros<mat>(nk,n);
  mat drsigma = zeros<mat>(nk,n);
  mat drLambda = zeros<mat>(nk,4);
  vec dre = zeros<vec>(nk);

  if(nprint>0) startMcmcTimer();

  for(int rep = 0; rep < ndpost; rep++) {
    
    cc = cgetC(e,k);
    bs = condd(Sigma);
    y = dy(y,x,cc,mu,as<mat>(bs["beta"]),as<vec>(bs["s"]),tau,sigma);
    
    //draw Sigma
    if(doSigma) {
      Res = y;
      Res.each_row() -= trans(mu);
      Res.each_col() -= tau;
      Res.each_col() /= sigma;
      
      S = trans(Res)*Res;
      rwout = rwishart(nu+n,solve(V+S,eye(p,p)));
      Sigma = as<mat>(rwout["IW"]);
    }
  
    //draw mu
    if(domu) {
      yd = y;
      yd.each_col() -= tau;
      Si = solve(Sigma,eye(p,p));
      Vmi = as_scalar(sum(1/pow(sigma,2)))*Si + Am;
      Rm = chol(Vmi);
      Ri = solve(trimatu(Rm),eye(p,p));
      Vm = solve(Vmi,eye(p,p));
      mm = Vm * (Si * (trans(yd) * (1/pow(sigma,2))) + Am * mubar);
      mu = vectorise(mm + Ri * as<vec>(rnorm(p)));
    }
      
    //draw tau
    if(dotau) {
      Ai = Lambda(0,0) - pow(Lambda(0,1),2)/Lambda(1,1);
      A = 1.0/Ai;
      onev = ones<mat>(p,1);
      Rm = chol(Sigma);
      xx = trans(solve(trans(Rm),onev));
      ytemp = trans(y);
      ytemp.each_col() -= mu;
      yy = trans(solve(trans(Rm),ytemp));
      xtx = accu(pow(xx,2)); //To get a sum of all the elements regardless of the argument type (ie. matrix or vector), use accu()
      xty = vectorise(xx*trans(yy));
      beta = A*Lambda(0,1)/Lambda(1,1);
      
      for(int j = 0; j<n; j++){
        s2 = xtx/pow(sigma[j],2) + A;
        s2 = 1.0/s2;
        m = s2*((xty[j]/pow(sigma[j],2)) + beta*(log(sigma[j])-Lambda(1,1)));
        tau[j] = m + sqrt(s2)*rnorm(1)[0];
      }
    }
     
    //draw sigma
    if(dosigma) {
      Rm = chol(Sigma);
      ytemp = y;
      ytemp.each_col() -= tau;
      ytemp = trans(ytemp);
      ytemp.each_col() -= mu;
      eps = solve(trans(Rm),ytemp);
      onesp = ones<rowvec>(p);
      ete = vectorise(onesp * pow(eps,2));
      
      a = Lambda(1,1);
      b = Lambda(0,1)/Lambda(0,0);
      s = sqrt(Lambda(1,1)-pow(Lambda(0,1),2)/Lambda(0,0));

      for(int j = 0; j<n; j++){
        pv = -(p+1)*log(gsigma) -.5*ete[j]/pow(gsigma,2) -.5*pow((log(gsigma)-(a+b*tau[j]))/s,2);
  	    pv = exp(pv-max(pv));
  	    pv = pv/sum(pv);
        //see http://gallery.rcpp.org/articles/using-the-Rcpp-based-sample-implementation/ for using sample
        sigma[j] = Rcpp::RcppArmadillo::sample(NumericVector(gsigma.begin(),gsigma.end()),1,false,NumericVector(pv.begin(),pv.end()))[0];
      }
    }
     
    //draw Lambda
    if(doLambda) {
      h = log(sigma);
      dat = join_rows(tau,h);        
      temp = cov(dat);
      moms << mean(tau) << mean(h) << temp(0,0) << temp(0,1) << temp(1,1); //element intialization
     
      SS = getS(Lambda,n,moms);
      rgl11 = gl11.elem(find(gl11 > pow(Lambda(0,1),2)/Lambda(1,1)));
      ng = rgl11.size();
      pv = zeros<vec>(ng);
      
      for(int j  = 0; j<ng; j++){
        Lambda(0,0) = rgl11[j];
        pv[j] = llL(Lambda,n,SS,VL,nuL);
      }
      
      pv = exp(pv-max(pv));
      pv = pv/sum(pv);
      Lambda(0,0) = Rcpp::RcppArmadillo::sample(NumericVector(rgl11.begin(),rgl11.end()),1,false,NumericVector(pv.begin(),pv.end()))[0];
      
      //cannot do multiple conditions per find() so it is done in two stages
      rgl12a = gl12.elem(find(gl12 < sqrt(Lambda(0,0)*Lambda(1,1))));
      rgl12 = rgl12a.elem(find(rgl12a > -sqrt(Lambda(0,0)*Lambda(1,1))));
      ng = rgl12.size();
      pv = zeros<vec>(ng);
      
      for(int j = 0; j<ng;j++){
        Lambda(0,1) = rgl12[j];
        Lambda(1,0) = Lambda(0,1);
        pv[j] = llL(Lambda,n,SS,VL,nuL);
      }
      
      pv = exp(pv-max(pv));
      pv = pv/sum(pv);
      
      Lambda(0,1) = Rcpp::RcppArmadillo::sample(NumericVector(rgl12.begin(),rgl12.end()),1,false,NumericVector(pv.begin(),pv.end()))[0];
      Lambda(1,0) = Lambda(0,1);
      
      rgl22 = gl22.elem(find(gl22 > pow(Lambda(0,1),2)/Lambda(0,0)));
      ng = rgl22.size();
      pv = zeros<vec>(ng);
      
      for(int j = 0;j<ng;j++){
        Lambda(1,1) = rgl22[j];
        SS = getS(Lambda,n,moms);
        pv[j] = llL(Lambda,n,SS,VL,nuL);
      }
      
      pv = exp(pv-max(pv));
      pv = pv/sum(pv);
      Lambda(1,1) = Rcpp::RcppArmadillo::sample(NumericVector(rgl22.begin(),rgl22.end()),1,false,NumericVector(pv.begin(),pv.end()))[0];
    }
    
    //draw e
    if(doe) {
      ng = ge.size();
      absege = abs(e-ge); 
      eiu = find(absege == min(absege));
      ei = eiu[0];
      
      if(ei == 1){
        pi = 2;
        qr = .5;
      } else if (ei == ng) {
        pi = ng-1;
        qr = .5;
      } else {
        pi = ei + rbinom(1,1,.5)[0]*2-1;
        qr = 1;
      }

      eold = ge[ei];
      eprop = ge[pi];
      
      llold = rlpx(x,eold,k,mu,tau,Sigma,sigma,ndghk);
      llprop = rlpx(x,eprop,k,mu,tau,Sigma,sigma,ndghk);
      lrat = llprop - llold + log(qr);
      
      if(lrat>0) {
        e = eprop;
      } else {
        minvec << 1 << exp(lrat);
        paccept = min(minvec);
        
        if(rbinom(1,1,paccept)[0]==1){
          e = eprop;
        } else {
          e = eold;
        }
      }
    }
     
    if (nprint>0) if((rep+1)%nprint==0) infoMcmcTimer(rep, R); 
  
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      drSigma(mkeep-1,span::all) = trans(vectorise(Sigma));
      drmu(mkeep-1,span::all) = trans(mu);
      drtau(mkeep-1,span::all) = trans(tau);
      drsigma(mkeep-1,span::all) = trans(sigma); 
      drLambda(mkeep-1,span::all) = trans(vectorise(Lambda));
      dre[mkeep-1] = e;
    }
  }
  
  if (nprint>0) endMcmcTimer();
  
  return List::create(
    Named("ndpost") = ndpost,
    Named("drmu") = drmu,
    Named("drtau") = drtau,
    Named("drsigma") = drsigma,
    Named("drLambda") = drLambda,
    Named("dre") = dre,
    Named("drSigma") = drSigma);
}
