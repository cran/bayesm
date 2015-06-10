#include "bayesm.h"
 
//FUNCTION SPECIFIC TO MAIN FUNCTION--------------------------------
vec rcomp(List comp) {
  
// Wayne Taylor 9/10/14

//purpose: draw multivariate normal with mean and variance given by comp 
// arguments:
//     comp is a list of length 2,
//     comp[[1]] is the mean and comp[[2]] is R^{-1} = comp[[2]], Sigma = t(R)%*%R

  vec mu = comp[0];
  mat rooti = comp[1];

  int dim = rooti.n_cols;
  mat root = solve(trimatu(rooti),eye(dim,dim)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  
  return(vectorise(mu+trans(root)*as<vec>(rnorm(mu.size()))));
}

//[[Rcpp::export]]
List rmixture(int n, vec pvec, List comps) {
                           
// Wayne Taylor 9/10/2014

// revision history:
//   commented by rossi 3/05
//
// purpose: iid draws from mixture of multivariate normals
// arguments:
//     n: number of draws
//     pvec: prior probabilities of normal components
//     comps: list, each member is a list comp with ith normal component
//                     ~N(comp[[1]],Sigma), Sigma = t(R)%*%R, R^{-1} = comp[[2]]
// output:
//  list of x (n by length(comp[[1]]) matrix of draws) and z latent indicators of
//  component

  //Draw vector of indices using base R 'sample' function
  mat prob(n,pvec.size());
  for(int i = 0; i<n; i++) prob(i,span::all) = trans(pvec);
  
  // Cumulatively add each row and take a uniform draw between 0 and the cumulative sum
  prob = cumsum(prob, 1);
  vec u = as<vec>(runif(n)) % prob.col(pvec.size()-1);
  
  // Evaluative each column of "prob" until the uniform draw is less than the cumulative value
  vec z = zeros<vec>(n);
  for(int i = 0; i<n; i++) while(u[i] > prob(i, z[i]++));

  List comp0 = comps[0];
  vec mu0 = comp0[0];
  mat x(n,mu0.size());
  
  //Draw from MVN from comp determined from z index
  //Note z starts at 1, not 0
  for(int i = 0; i<n; i++) {
    x(i,span::all) = trans(rcomp(comps[z[i]-1]));
  }

  return List::create(
    Named("x") = x, 
    Named("z") = z);
}
