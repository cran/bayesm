#include "bayesm.h"
 
// [[Rcpp::export]]
vec cgetC(double e, int k){

//Wayne Taylor 4/29/15

//purpose: get a list of cutoffs for use with scale usage problems

//arguments:
//   e: the "e" parameter from the paper
//   k: the point scale, eg. items are rated from 1,2,...k
// output:
//   vector of grid points

  vec temp = zeros<vec>(k-1);
  for(int i = 0; i<(k-1); i++) temp[i] = i + 1.5;
  double m1 = sum(temp);
  temp = pow(temp,2);
  double m2 = sum(temp);
    
  vec c = zeros<vec>(k+1);
  
  //first sum to get s's, this is a waste since it should be done
  //once but I don't want to see this things anywhere else and it should take no time
  double s0 = k-1;
  double s1=0.0,s2=0.0,s3=0.0,s4=0.0;
  for(int i=1;i<k;i++) {s1+=i; s2+=i*i; s3+= i*i*i; s4+=i*i*i*i;}

  // now make quadratic for b (just as in Peter's code)
  double aq = s0*s2-s1*s1;
  double bq = 2*e*s0*s3-2*e*s1*s2;
  double cq = m1*m1 - m2*s0 + e*e*s0*s4 - e*e*s2*s2;

  //get a and b
  double det = bq*bq - 4*aq*cq;
  if(det<0) stop("no solution for c's given e and m1, m2 \n");
  double b=(-bq+sqrt(det))/(2.0*aq);
  double a=(m1-b*s1-e*s2)/s0;

  //make c
  c[0]= -1000.0;
  c[k]= 1000.0;
  for(int i=1;i<k;i++) c[i] = a+b*i+e*i*i;
   
  return(sort(c));
}
