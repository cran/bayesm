rmultiregfp=
function(Y,X,Fparm)
{
#
# revision history:
#    changed 1/11/05 by P. Rossi to fix sum of squares error
#
# purpose:
#    draw from posterior for Multivariate Regression Model
# arguments:
#    Y is n x m matrix
#    X is n x k
#    Fparm is a list of fixed parameters created by init.rmultiregfp    
# output:
#    list of B, Sigma draws of matrix of coefficients and Sigma matrix
# model:
#    Y=XB+U  cov(u_i) = Sigma
#    B is k x m matrix of coefficients
# priors:  beta|Sigma  ~ N(betabar,Sigma (x) A^-1)
#                   betabar=vec(Bbar)
#                   beta = vec(B) 
#          Sigma ~ IW(nu,V) or Sigma^-1 ~ W(nu, V^-1)
n=nrow(Y)
m=ncol(Y)
k=ncol(X)
#
# first draw Sigma
#
W=rbind(X,Fparm$RA)
Z=rbind(Y,Fparm$RABbar)
#   note:  Y,X,A,Bbar must be matrices!
IR=Fparm$IR
#                      W'W = R'R  &  (W'W)^-1 = IRIR'  -- this is the UL decomp!
Btilde=crossprod(t(IR))%*%crossprod(W,Z)   
#                      IRIR'(W'Z) = (X'X+A)^-1(X'Y + ABbar)
S=crossprod(Z-W%*%Btilde)
#                      E'E
rwout=rwishart(Fparm$nu+n,chol2inv(chol(Fparm$V+S)))
#
# now draw B given Sigma
#   note beta ~ N(vec(Btilde),Sigma (x) Covxxa)
#       Cov=(X'X + A)^-1  = IR t(IR)  
#       Sigma=CICI'    
#       therefore, cov(beta)= Omega = CICI' (x) IR IR' = (CI (x) IR) (CI (x) IR)'
#	so to draw beta we do beta= vec(Btilde) +(CI (x) IR)vec(Z_mk)  
#			Z_mk is m x k matrix of N(0,1)
#	since vec(ABC) = (C' (x) A)vec(B), we have 
#		B = Btilde + IR Z_mk CI'
#
B = Btilde + IR%*%matrix(rnorm(m*k),ncol=m)%*%t(rwout$CI)
return(list(B=B,Sigma=rwout$IW))
}

