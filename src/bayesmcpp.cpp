// R.McCulloch, 12/04  code for scale usage R function (rScaleUsage)
#include <iostream>
#include <algorithm>

extern "C" {
#include <R.h>
#include <Rmath.h>
};

extern "C" void getC(double *ep,int *kp, double *m1p, double *m2p,  double *c);
extern "C" void dy(int *p, int *nob,  double *y, int *x, double *c, double *mu, double *beta, double *s, double *tau, double *sigma);

void getC(double *ep,int *kp, double *m1p, double *m2p, double *c)
{
   double e = *ep;
   int k = *kp;
   double m1 = *m1p;
   double m2 = *m2p;

   //first sum to get s's, this is a waste since it should be done
   //once but I don't want to see this things anywhere else and it should take no time
   double s0 = (double)(k-1);
   double s1=0.0,s2=0.0,s3=0.0,s4=0.0;
   for(int i=1;i<k;i++) {s1+=i; s2+=i*i; s3+= i*i*i; s4+=i*i*i*i;}

   // now make quadratic for b (just as in Peter's code)
   double aq = s0*s2-s1*s1;
   double bq = 2*e*s0*s3-2*e*s1*s2;
   double cq = m1*m1 - m2*s0 + e*e*s0*s4 - e*e*s2*s2;

   //get a and b
   double det = bq*bq - 4*aq*cq;
   if(det<0) std::cout << "error: no solution for c's given e and m1, m2" << std::endl;
   double b=(-bq+sqrt(det))/(2.0*aq);
   double a=(m1-b*s1-e*s2)/s0;

   //make c
   c[0]= -1000.0;
   c[k]= 1000.0;
   for(int i=1;i<k;i++) c[i] = a+b*i+e*i*i;
   
   std::sort(c,c+k+1);

}



void d1y(int p, double *y, int *x, double *c, double *mu, double *beta, double *s, double tau, double sigma)
{
   //std::cout << "int main of d1y" << std::endl;

   GetRNGstate();
   double cm,cs; //cm = conditional mean, cs = condtional standard deviation
   double u;    // uniform for truncated normal draw
   double a,b;  // standardized truncation points
   double pa,pb; // cdf at truncation points

   //loop over coordinates of y
   for(int i=0;i<p;i++) {
      //compute conditonal mean and standard deviation
      cs = s[i]*sigma;
      cm = mu[i]+tau;
      for(int j=0;j<i;j++) cm += (*(beta+i*(p-1)+j))*(y[j]-mu[j]-tau);
      for(int j=(i+1);j<p;j++) cm += (*(beta+i*(p-1)+j-1))*(y[j]-mu[j]-tau);
      //draw truncated normal
      // y~N(cm,cs^2) I[c[x[i]-1],c[x[i])
      a = (c[x[i]-1]-cm)/cs;  b = (c[x[i]]-cm)/cs;
      pa = pnorm(a,0.0,1.0,1,0); pb = pnorm(b,0.0,1.0,1,0);
      u = unif_rand();
      y[i] = cm + cs*qnorm(u*pb + (1-u)*pa,0.0,1.0,1,0);
   }
   PutRNGstate();
}

void dy(int *p, int *nob, double *y, int *x, double *c, double *mu, double *beta, double *s, double *tau, double *sigma)
{
   for(int i=0;i<(*nob);i++) {
      d1y(*p,y+i*(*p),x+i*(*p),c,mu,beta,s,*(tau+i),*(sigma+i));
   }
}
