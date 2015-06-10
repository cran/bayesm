#include "bayesm.h"
 
//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------
mat ztoSim(vec const& z){

// function to convert indicator vector to Similarity matrix
// Sim is n x n matrix, Sim[i,j]=1 if pair(i,j) are in same group
// z is n x 1 vector of indicators (1,...,p)

  int n = z.size();
  vec onevec = ones<vec>(n); // equivalent to zvec=c(rep(z,n)) in R
  vec zvec = kron(onevec, z);
  vec zcomp = kron(z, onevec);// equivalent to as.numeric((zvec==zcomp)) in R
  mat Sim = zeros<mat>(n*n,1);
  
  for (int i=0; i<n*n; i++){
    if (zvec[i]==zcomp[i]) Sim(i,0) = 1;
  }
  
  Sim.reshape(n,n);
  
  return (Sim);
}

vec Simtoz(mat const& Sim){

// function to convert Similarity matrix to indicator vector
// Sim is n x n matrix, Sim[i,j]=1 if pair(i,j) are in same group
// z is vector of indicators from (1,...,p) of group memberships (dim n)

  int count, i, j;
  int n = Sim.n_cols;
  vec z = zeros<vec>(n);  
  int groupn = 1;
  
  for (i=0; i<n; i++){
    count = 0;
    for (j=0; j<n; j++){    
      if ((z[j]==0) & (Sim(j,i)==1)){
        z[j] = groupn;
        count = count + 1;
      }      
    }
    if (count>0){
      groupn = groupn + 1;
    }
  }
  
  return (z);
} 

//MAIN FUNCTION---------------------------------------------------------------------------------------
// [[Rcpp::export]]
List clusterMix_rcpp_loop(mat const& zdraw, double cutoff, bool SILENT, int nprint){

// Keunwoo Kim 10/06/2014

// Purpose: 
//    cluster observations based on draws of indicators of
//    normal mixture components

// Arguments:
//    zdraw is a R x nobs matrix of draws of indicators (typically output from rnmixGibbs)
//    the rth row of zdraw contains rth draw of indicators for each observations
//    each element of zdraw takes on up to p values for up to p groups. The maximum
//    number of groups is nobs.  Typically, however, the number of groups will be small
//    and equal to the number of components used in the normal mixture fit.

//    cutoff is a cutoff used in determining one clustering scheme it must be 
//    a number between .5 and 1.

//    nprint - print every nprint'th draw

// Output: 
//    two clustering schemes each with a vector of length nobs which gives the assignment
//    of each observation to a cluster

//    clustera (finds zdraw with similarity matrix closest to posterior mean of similarity)
//    clusterb (finds clustering scheme by assigning ones if posterior mean of similarity matrix cutoff and computing associated z )
 
  int rep, i;
  uword index; // type uword means unsigned integer. Necessary for finding the index of min.    
  int nobs = zdraw.n_cols;    
  char buf[32];
  
  // compute posterior mean of Similarity matrix
  if (!SILENT){
    Rcout << "Computing Posterior Expectation of Similarity Matrix\n";
    Rcout << "processing draws ...\n";
  }
  
  mat Pmean = zeros<mat>(nobs, nobs);
  int R = zdraw.n_rows;
  
  for (rep=0; rep<R; rep++){
    Pmean = Pmean + ztoSim(trans(zdraw(rep,span::all)));
    if (!SILENT){
      if ((rep+1)%nprint==0){        
        sprintf(buf, "  %d\n", rep+1);
        Rcout <<  buf;
      }
    }
  }
  
  Pmean = Pmean/R;
  
  // now find index for draw which minimizes discrepancy between
  // post exp of similarity and sim implied by that z
  if (!SILENT){
    Rcout << " \n";
    Rcout << "Look for zdraw which minimizes loss \n";
    Rcout << "processing draws ... \n";
  }
  
  vec loss = zeros<vec>(R);
  
  for (rep=0; rep<R; rep++){
    loss[rep] = accu(abs(Pmean-ztoSim(trans(zdraw(rep,span::all)))));
    if (!SILENT){
      if ((rep+1)%nprint==0){
        sprintf(buf, "  %d\n", rep+1);
        Rcout <<  buf;
      }
    }
  }
  
  loss.min(index);  
  vec clustera = trans(zdraw(index,span::all));
  
  // now due clustering by assigning Similarity to any (i,j) pair for which
  // Pmean > cutoff
  vec Pmeanvec = vectorise(Pmean);
  mat Sim = zeros<mat>(nobs*nobs,1);
  
  for (i=0; i<nobs*nobs; i++){
    if (Pmeanvec[i]>=cutoff) Sim(i,0) = 1;      
  }
  
  Sim.reshape(nobs,nobs);  
  vec clusterb = Simtoz(Sim);

  return List::create(
      Named("clustera") = clustera,
      Named("clusterb") = clusterb);
}
