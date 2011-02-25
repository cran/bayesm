#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <Rdefines.h>

/* modified by rossi 7/06 to remove thetaStar argument and copy */
/* modified by rossi 7/06 to remove theta copy and modify directly */
/* modified by rossi 8/10 to fix error in allocation of size of newrow */

/* function to make multinomial draws */

int rmultin(double *probs, int nprob )
{
	double cumprob,rnd;
	int i;
	GetRNGstate();
	rnd=unif_rand();
	cumprob=0.0;
	for (i=0; i < nprob; i++){
		if(rnd > cumprob && rnd <= cumprob+probs[i])
		{ break; }
		else
		{cumprob=cumprob+probs[i];}
	}
	PutRNGstate(); 
	return i+1 ;
}
		

/* gets a row of a matrix and returns this as a matrix by setting dim */

SEXP getrow(SEXP mat, int row, int nrow, int ncol){
   int i,ind;
   SEXP ans, ndim;
   PROTECT(ans=NEW_NUMERIC(ncol));
   PROTECT(ndim=NEW_INTEGER(2));
   for(i =0; i < ncol; i++){
	   ind=i*nrow+row;
	   NUMERIC_POINTER(ans)[i]=NUMERIC_POINTER(mat)[ind];
   }
   INTEGER_POINTER(ndim)[0]=1;
   INTEGER_POINTER(ndim)[1]=ncol;
   SET_DIM(ans,ndim);
   UNPROTECT(2);
return(ans);
}


/* theta draw routine to be used with .Call */

SEXP  thetadraw( SEXP y,  SEXP ydenmatO, SEXP indicO, SEXP q0v, SEXP p,
	SEXP theta,  SEXP lambda, SEXP eta,
                  SEXP thetaD, SEXP yden,
		  SEXP maxuniqS,SEXP nuniqueS,
                  SEXP rho) {
   int nunique,n,ncol,j,i,maxuniq,inc,index,ii,jj,ind ;
   SEXP R_fc_thetaD, R_fc_yden, yrow, ydim, onetheta, lofone, newrow,
	ydenmat, ydendim ;
   double *probs;
   int *indmi;
   int *indic;
   double sprob;

   nunique=INTEGER_VALUE(nuniqueS);
   n=length(theta);
   maxuniq=INTEGER_VALUE(maxuniqS);

   /* create new lists for use and output */ 
   PROTECT(lofone=NEW_LIST(1));
  
   /* create R function call object, lang4 creates a pairwise (linked) list with
      4 values -- function, first arg, sec arg, third arg.  R_NilValue is a placeholder until
      we associate first argument (which varies in our case) */
   PROTECT(R_fc_thetaD=lang4(thetaD,R_NilValue,lambda,eta));
   PROTECT(R_fc_yden=lang4(yden,R_NilValue,y,eta));

   PROTECT(ydim=GET_DIM(y));
   ncol=INTEGER_POINTER(ydim)[1];
   PROTECT(yrow=NEW_NUMERIC(ncol));
   PROTECT(newrow=NEW_NUMERIC(n));
   PROTECT(ydenmat=NEW_NUMERIC(maxuniq*n));
   PROTECT(ydendim=NEW_INTEGER(2));
   INTEGER_POINTER(ydendim)[0]=maxuniq;
   INTEGER_POINTER(ydendim)[1]=n;

   /* copy iformation from R objects that will be modified      
      note that we must access elements in the lists (generic vectors) by using VECTOR_ELT
      we can't use the pointer and deferencing directly like we can for numeric and integer
      vectors */
   for(j=0;j < maxuniq*n; j++){NUMERIC_POINTER(ydenmat)[j]=NUMERIC_POINTER(ydenmatO)[j];}
   SET_DIM(ydenmat,ydendim); 

   /* allocate space for local vectors */
   probs=(double *)R_alloc(n,sizeof(double));
   indmi=(int *)R_alloc((n-1),sizeof(int));
   indic=(int *)R_alloc(n,sizeof(int));
   
   /* copy information from R object indicO to indic */
   for(j=0;j < n; j++) {indic[j]=NUMERIC_POINTER(indicO)[j];}

   /* start loop over observations */

   for(i=0;i < n; i++){
	 probs[n-1]=NUMERIC_POINTER(q0v)[i]*NUMERIC_POINTER(p)[n-1];

	 /* make up indmi -- vector of length n-1 consisting of -i as in R notation --
	    1, ...,i-1, ,i+1,...,n */
	 inc=0;
	 for(j=0;j < (n-1); j++){
		 if(j==i) {inc=inc+1;};
		 indmi[j]=inc;
		 inc=inc+1;
	 }
	 for(j=0;j < (n-1); j++){
		 ii=indic[indmi[j]]; jj=i;      /* find element ydenmat(ii,jj+1) */
		 index=jj*maxuniq+(ii-1);
		 probs[j]=NUMERIC_POINTER(p)[j]*NUMERIC_POINTER(ydenmat)[index];
	 }
	 sprob=0.0;
	 for(j=0;j<n;j++){sprob=sprob+probs[j];}
	 for(j=0;j<n;j++){probs[j]=probs[j]/sprob;}
	 ind=rmultin(probs,n);
          
	 if(ind == n){
                 yrow=getrow(y,i,n,ncol);
                 SETCADR(R_fc_thetaD,yrow);   /* set the second argument to yrow -- head of the tail */
        	 onetheta=eval(R_fc_thetaD,rho);
                 SET_ELEMENT(theta,i,onetheta);
		 if((nunique) > (maxuniq-1)) {error("max number of unique thetas exceeded");}
		                             /* check to make sure we don't exceed max number of unique theta */
	         SET_ELEMENT(lofone,0,onetheta);
	         SETCADR(R_fc_yden,lofone);
	         newrow=eval(R_fc_yden,rho);
	         for(j=0;j<n; j++)
		    { NUMERIC_POINTER(ydenmat)[j*maxuniq+nunique]=NUMERIC_POINTER(newrow)[j];}
		 indic[i]=nunique+1;
		 nunique=nunique+1;
	 }
	 else {
		 onetheta=VECTOR_ELT(theta,indmi[ind-1]);
		 SET_ELEMENT(theta,i,onetheta);
		 indic[i]=indic[indmi[ind-1]];
	 }
    }

    UNPROTECT(8);
    return(nuniqueS);     /* returns argument -- function now is called for its effect on theta */
}
 

