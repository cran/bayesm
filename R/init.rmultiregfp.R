init.rmultiregfp=
function(X,A,Bbar,nu,V) 
{
#
#       revision history:
#       revised 1/11/05 by P. Rossi
# purpose: 
#    prepare parameters for rmultiregfp  
#        rmultiregfp is designed to take advantage of the fact that posterior calculations
#        don't change X or prior parms even when the "data" or Y is changing
#
#        use to create list for call of rmultiregfp as in:
#               Fparm=init.rmultiregfp(X,A,Bbar,nu,V)
#	        rmultiregfp(Y,X,Fparm)
#
# 
# for multivariate regression model Y = XB + U  Y is n x m, X, n x k, B is k x m
#			ith row of U is N(0,Sigma)
#
#	 prior parms:  vec(B)| Sigma ~ N(vec(Bbar),Sigma (x) A^-1)  
#                      Sigma ~ IW(nu,V)
#       
k=ncol(X)
IR=backsolve(chol(crossprod(X)+A),diag(k))
RA=chol(A)
RABbar=RA%*%Bbar
list(IR=IR,RA=RA,RABbar=RABbar,nu=nu,V=V)
}
