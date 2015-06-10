#include "bayesm.h"
 
//W. Taylor: we considered moving the output to struct formats but the efficiency
//  gains were limited and the conversions back and forth between Lists and struct were cumbersome

List drawCompsFromLabels(mat const& y,  mat const& Bbar, 
                         mat const& A, int nu, 
                         mat const& V,  int ncomp,
                         vec const& z){
                           
// Wayne Taylor 3/18/2015

// Function to draw the components based on the z labels
  
  vec b, r, mu;
  mat yk, Xk, Ck, sigma, rooti, S, IW, CI;
  List temp, rw, comps(ncomp);
  
  int n = z.n_rows;
  vec nobincomp = zeros<vec>(ncomp);
  
  //Determine the number of observations in each component
  for(int i = 0; i<n; i++) {
    nobincomp[z[i]-1]++; //Note z starts at 1, not 0
  }
  
  //Draw comps
  for(int k = 0; k<ncomp; k++){
    
    if(nobincomp[k] > 0) {
      // If there are observations in this component, draw from the posterior
      
      yk = y.rows(find(z==(k+1))); //Note k starts at 0 and z starts at 1
      Xk = ones(nobincomp[k], 1);

      temp = rmultireg(yk, Xk, Bbar, A, nu, V);
      
      sigma = as<mat>(temp["Sigma"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      rooti = solve(trimatu(chol(sigma)),eye(sigma.n_rows,sigma.n_cols)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
      
      mu = as<vec>(temp["B"]);

      comps(k) = List::create(
        Named("mu") = NumericVector(mu.begin(),mu.end()), //converts to a NumericVector, otherwise it will be interpretted as a matrix
        Named("rooti") = rooti
      );
      
    } else {
      // If there are no obervations in this component, draw from the prior
      S = solve(trimatu(chol(V)),eye(V.n_rows,V.n_cols));
      S = S * trans(S); 
      
      rw = rwishart(nu, S);
      
      IW = as<mat>(rw["IW"]);
      CI = as<mat>(rw["CI"]);
      
      rooti = solve(trimatu(chol(IW)),eye(IW.n_rows,IW.n_cols));        
      b = vectorise(Bbar);
      r = rnorm(b.n_rows,0,1);
      
      mu = b + (CI * r) / sqrt(A(0,0));
  	
		  comps(k) = List::create(
			  Named("mu") = NumericVector(mu.begin(),mu.end()), //converts to a NumericVector, otherwise it will be interpretted as a matrix
			  Named("rooti") = rooti);
    } 
  }

  return(comps);
}

vec drawLabelsFromComps(mat const& y, vec const& p, List comps) {
  
// Wayne Taylor 3/18/2015

// Function to determine which label is associated with each y value
  
  double logprod;
  vec mu, u;
  mat rooti;
  List compsk;
  
  int n = y.n_rows;
  vec res = zeros<vec>(n);
  int ncomp  = comps.size();
  mat prob(n,ncomp);
  
  for(int k = 0; k<ncomp; k++) {
    compsk = comps[k];
    mu = as<vec>(compsk["mu"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    rooti = as<mat>(compsk["rooti"]);

    //Find log of MVN density using matrices
    logprod = log(prod(diagvec(rooti)));
    mat z(y);
    z.each_row() -= trans(mu); //subtracts mu from each row in z
    z = trans(rooti) * trans(z);
    z = -(y.n_cols/2.0) * log(2*M_PI) + logprod - .5 * sum(z % z, 0); // operator % performs element-wise multiplication
      
    prob.col(k) =  trans(z);
  }

  prob = exp(prob);
  prob.each_row() %= trans(p); //element-wise multiplication

  // Cumulatively add each row and take a uniform draw between 0 and the cumulative sum
  prob = cumsum(prob, 1);
  u = as<vec>(runif(n)) % prob.col(ncomp-1);
  
  // Evaluative each column of "prob" until the uniform draw is less than the cumulative value
  for(int i = 0; i<n; i++) {
    while(u[i] > prob(i, res[i]++));
  }
  
  return(res);
}

vec drawPFromLabels(vec const& a, vec const& z) {
  
// Wayne Taylor 9/10/2014

// Function to draw the probabilities based on the label proportions
  
  vec a2 = a;
  int n = z.n_rows;
  
  //Count number of observations in each component
  for(int i = 0; i<n; i++) a2[z[i]-1]++; //z starts at 1, not 0
  
  return rdirichlet(a2);
}

//[[Rcpp::export]]
List rmixGibbs( mat const& y,  mat const& Bbar, 
                mat const& A, int nu, 
                mat const& V,  vec const& a, 
                vec const& p,  vec const& z) {

// Wayne Taylor 9/10/2014

/*
    // Revision History: R. McCulloch 11/04 P. Rossi 3/05 put in
    // backsolve and improved documentation
    // 
    // purpose: do gibbs sampling inference for a mixture of
    // multivariate normals
    // 
    // arguments: y: data, rows are observations, assumed to be iid
    // draws from normal mixture Bbar,A,nu,V: common prior for mean
    // and variance of each normal component
    // 
    // note: Bbar should be a matrix. usually with only one row
    // 
    // beta ~ N(betabar,Sigma (x) A^-1) betabar=vec(Bbar) Sigma ~
    // IW(nu,V) or Sigma^-1 ~ W(nu, V^-1) note: if you want Sigma ~
    // A, use nu big and rwishart(nu,nu(A)^{-1})$IW a: Dirichlet
    // parameters for prior on p p: prior probabilities of normal
    // components z: components indentities for each observation
    // (vector of intergers each in {1,2,...number of components})
    // comps: list, each member is a list comp with ith normal
    // component ~N(comp[[1]],Sigma), Sigma = t(R)%*%R, R^{-1} =
    // comp[[2]] Output: list with elements [[1]=$p, [[2]]=$z, and
    // [[3]]=$comps, with the updated values

    */
  
  List comps = drawCompsFromLabels(y, Bbar, A, nu, V, a.size(), z);
  
  vec z2 = drawLabelsFromComps(y, p, comps);
  
  vec p2 = drawPFromLabels(a, z2);

  return List::create(
    Named("p") = p2,
    Named("z") = z2,
    Named("comps") = comps);
}
