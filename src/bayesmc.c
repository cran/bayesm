#include <R.h>
#include <Rmath.h>
#include <math.h>

void condmom(double *x, double *mu, double *sigi, int p, int j, double *m, double *csig)
{
/*	function to compute moments of x[j] | x[-j]  */

	int ind,i,jm1;
	double csigsq;
	jm1=j-1;
	ind = p*jm1;
	csigsq = 1./sigi[ind+jm1];

	*m = 0.0;
	for (i=0 ; i < p; ++i)
		{
		if (i != jm1) 
			{*m +=  - csigsq*sigi[ind+i]*(x[i]-mu[i]);}
		}
	*m=mu[jm1]+*m ;
	*csig=sqrt(csigsq);
}

double rtrun(double mu, double sigma,double trunpt, int above) 
{
/*	function to draw truncated normal
		above=1 means from above b=trunpt, a=-inf
		above=0 means from below a=trunpt, b= +inf   
		modified by rossi 6/05 to check arg to qnorm
*/
	double FA,FB,rnd,result,arg ;
	if (above) {
		FA=0.0; FB=pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
			}
	else {
		FB=1.0; FA=pnorm(((trunpt-mu)/(sigma)),0.0,1.0,1,0);
		}
	
	GetRNGstate();
	rnd=unif_rand();
	arg=rnd*(FB-FA)+FA;
	if(arg > .999999999) arg=.999999999;
	if(arg < .0000000001) arg=.0000000001;
	result = mu + sigma*qnorm(arg,0.0,1.0,1,0);
	PutRNGstate();
	return result;
}


void drawwi(double *w, double *mu, double *sigmai,int *p, int *y)
{
/*	function to draw w_i by Gibbing's thru p vector   */

	int i,j,above;
	double bound;
	double mean, csig;

		for (i=0; i < *p; ++i) 
		{	
			bound=0.0;
		    	for (j=0; j < *p ; ++j) 
				{ if (j != i) {bound=fmax2(bound,w[j]); }}
			if (*y == i+1) 	
				above = 0;
			else 
				above = 1;

		condmom(w,mu,sigmai,*p,(i+1),&mean,&csig);
		w[i]=rtrun(mean,csig,bound,above);

		}
}

void draww(double *w, double *mu, double *sigmai, int *n, int *p, int *y) 
{
/*	function to gibbs down entire w vector for all n obs  */
	int i, ind;
	for (i=0; i < *n; ++i)
	{
		ind= *p * i;
		drawwi(w+ind,mu+ind,sigmai,p,y+i);
	}
}


void drawwi_mvp(double *w, double *mu, double *sigmai,int *p, int *y)
{
/*	function to draw w_i for Multivariate Probit  */

	int i,above;
	double mean, csig;

		for (i=0; i < *p; ++i) 
		{	
			if (y[i]) 	
				above = 0;
			else 
				above = 1;

		condmom(w,mu,sigmai,*p,(i+1),&mean,&csig);
		w[i]=rtrun(mean,csig,0.0,above);

		}
}

void draww_mvp(double *w, double *mu, double *sigmai, int *n, int *p, int *y) 
{
/*	function to gibbs down entire w vector for all n obs  */
	int i, ind;
	for (i=0; i < *n; ++i)
	{
		ind= *p * i;
		drawwi_mvp(w+ind,mu+ind,sigmai,p,y+ind);
	}
}

double root(double c1, double  c2, double *tol,int *iterlim)
{
/*	function to find root of c1 - c2u = lnu */
   int iter;
   double uold, unew;
   uold=1.;
   unew=0.00001;
   iter=0;
   while (iter <= *iterlim && fabs(uold-unew) > *tol )
      {
      uold=unew;
      unew=uold + (uold*(c1 -c2*uold -  log(uold)))/(1. + c2*uold); 
      if(unew < 1.0e-50) unew=1.0e-50;
      iter=iter+1;
      }
   return unew;
}
   
void callroot(int *n,double *c1, double *c2, double *tol, int *iterlim,double *u)
{
   int i;
   for (i=0;i < *n; ++i)
   {
	u[i]=root(c1[i],c2[i],tol,iterlim);
   }
}



void ghk_oneside(double *L, double* trunpt, int *above, int *dim, int *n, double *res)
/*	routine to implement ghk with a region defined by truncation only on one- side
 						r mcculloch 8/04
        if above=1, then we truncate component i from above at point trunpt[i-1]
        L is lower triangular root of Sigma
	random vector is assumed to have zero mean
    	n is number of draws to use in GHK	
	modified 6/05 by rossi to check arg into qnorm
*/
{
   int i,j,k;
   double mu,tpz,u,prod,pa,pb,arg;
   double *z;
   z = (double *)R_alloc(*dim,sizeof(double));
   GetRNGstate();
   *res = 0.0;
   for(i=0;i<*n;i++) {
      prod=1.0;
      for(j=0;j<*dim;j++) {
         mu=0.0; for(k=0;k<j;k++) mu += L[k*(*dim)+j]*z[k];
	 tpz = (trunpt[j]-mu)/L[j*(*dim)+j];
	 if(above[j]) {
	    pa=0.0; pb = pnorm(tpz,0.0,1.0,1,0);   
	 }
	 else {
	    pb=1.0; pa = pnorm(tpz,0.0,1.0,1,0);
	 }
	 prod *= pb-pa;
	 u = unif_rand();
	 arg=u*pb+(1.-u)*pa;
	 if(arg > .999999999) arg=.999999999;
	 if(arg < .0000000001) arg=.0000000001;
	 z[j] = qnorm(arg,0.0,1.0,1,0);
      }
      *res += prod;
   }
   *res /= (double)(*n);
   PutRNGstate();
}
void ghk(double *L, double* a, double *b, int *dim, int *n, double *res)
/*	routine to implement ghk with a region : a[i-1] <= x_i <= b[i-1]
 						r mcculloch 8/04
        L is lower triangular root of Sigma
	random vector is assumed to have zero mean
    	n is number of draws to use in GHK	
	modified 6/05 by rossi to check arg into qnorm
*/
{
   int i,j,k;
   double aa,bb,pa,pb,u,prod,mu,arg;
   double *z;
   z = (double *)R_alloc(*dim,sizeof(double));
   GetRNGstate();
   *res=0.0;
   for(i=0;i<*n;i++) {
      prod = 1.0;
      for(j=0;j<*dim;j++) {
         mu=0.0; for(k=0;k<j;k++) mu += L[k*(*dim)+j]*z[k];
	 aa=(a[j]-mu)/L[j*(*dim)+j]; bb = (b[j]-mu)/L[j*(*dim)+j];
	 pa = pnorm(aa,0.0,1.0,1,0); pb = pnorm(bb,0.0,1.0,1,0);
	 prod *= pb-pa;
	 u = unif_rand();
	 arg=u*pb+(1.-u)*pa;
	 if(arg > .999999999) arg=.999999999;
	 if(arg < .0000000001) arg=.0000000001;
	 z[j] = qnorm(arg,0.0,1.0,1,0);
      }
      *res += prod;
   }
   *res /= (double)(*n);
   PutRNGstate();
}

void ghk_vec(int *n,double *L, double *trunpt,int *above, int *dim, int *r, double *res)
{
/* routine to call ghk_oneside for n different truncation points stacked in to the
   vector trunpt  -- puts n results in vector res
                                          p rossi 12/04
*/
	int i, ind;
	for (i=0; i < *n; ++i)
	{
		ind = *dim * i;
		ghk_oneside(L,trunpt + ind,above,dim,r,res+i);
	}
}

void cuttov(double *ut,double *v, int *dim)
/*
purpose: write upper triangular (ut) to vector (v), goes down columns, omitting zeros
arguments:
   ut: upper triangular matrix, stored as series of columns (including the zeros)
   v: vector ut is copied to, on input must have correct length
   dim: ut is dim x dim, v is dim*(dim+1)/2
*/
{
   int ind=0;
   int i,j;
   for(i=0;i<(*dim);i++) {
      for(j=0;j<=i;j++) {
         v[ind] = ut[i*(*dim)+j];
         ind += 1;
      }
   }
}
void cvtout(double *v, double *ut, int *dim)
/*
purpose: write vector (v) to upper triangular (inverse of cuttov above)
arguments:
   v: vector
   ut: upper triangulare matrix, columns stacked, zeros included
   dim: ut is dim x dim, v is dim*(dim+1)/2
*/
{
   int ind=0;
   int i,j;
   for(i=0;i<(*dim);i++) {
      for(j=(i+1);j<(*dim);j++) ut[i*(*dim)+j]=0.0;
      for(j=0;j<=i;j++) {
         ut[i*(*dim)+j] = v[ind];
         ind += 1;
      }
   }
}
void clmvn(double *x, double *mu, double *riv, int *dim, double *res)
/*
purpose:
   calculate log of multivariate density 
   evaluated at x
   mean is mu, and covariance matrix t(R)%*%R and riv is vector version of the inverse of R
arguments:
   x: compute log(f(x))
   mu, riv: x~N(mu,t(R)%*%R), riv is vector version of R^{-1}
   dim: dimension of x
   res: place to put result
*/
{
   int i,j;
   double sum = 0.0;
   double prod = 1.0;
   double z;
   int ind = 0;
   for(i=0;i<(*dim);i++) {
      z = 0.0;
      for(j=0;j<=i;j++) {z += riv[ind]*(x[j]-mu[j]); ind += 1;}
      sum += z*z;
      prod *= riv[ind-1];
   }
   *res = log(prod) -.5*sum;
}
void crdisc(double *p, int *res)
/*
purpose: draw from a discrete distribution
arguments:
   p: vector of probabilities
   res: draw is in {1,2,...length(p)}, giving the draw's category
*/
{
   double u,sum;
   GetRNGstate();
   u = unif_rand();
   *res = 1;
   sum = p[*res -1];
   while(sum<u) {sum += p[*res];(*res) +=1;}
   PutRNGstate();
}
void crcomp(double *x, double *mu, double *riv, double *p, int *dim, int *nc, int *res)
/*
purpose: draw component of x, where x is drawn form mixture of multivariate normal components
arguments:
   x: observed vector x
   mu: matrix of class means, each column gives a mean vector
   riv: matrix of class covariances, each column gives a vector version of R^{-1}, Sigma = t(R)%*%R
   p: prior class probabilities
   dim: dimension of x (and mu, and Sigma)
   nc: number of classes
   res: result
note:
   mu is column stacked version of a dim x nc matrix
   riv is column stacked version of a dim*(dim+1)/2 x nc matrix
*/
{
   double *post;
   double max,sum;
   int dim_riv = (*dim)*((*dim)+1)/2;
   int i;
   post = (double *)R_alloc(*nc,sizeof(double));
   clmvn(x,mu,riv,dim,post);
   max = *post;
   for(i=1;i<(*nc);i++) {
      clmvn(x,mu+i*(*dim),riv+i*dim_riv,dim,post+i);
      if(*(post+i) > max) max = *(post+i);
   }
   sum = 0.0;
   for(i=0;i<(*nc);i++) { post[i] = exp(post[i]-max)*p[i]; sum += post[i];}
   for(i=0;i<(*nc);i++) post[i] /= sum;
   crdisc(post,res);
}
void crcomps(double *x, double *mu, double *riv, double *p, int *dim, int *nc, int *nob, int *res)
/*
purpose: x represents a matrix, whose columns are draws from a normal mixture, draw component membership for each x
arguments:
   all the same as crcomp, except x is now column stacked version of dim x nob matrix
   and nob is the number of observations
   res is now of length nob
*/
{
   int i;
   for(i=0;i<(*nob);i++) {
      crcomp(x+i*(*dim),mu,riv,p,dim,nc,res+i);
   }
}

