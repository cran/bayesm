rscaleUsage=
function(Data,Prior,Mcmc) 
{
#
# purpose: run scale-usage mcmc
#    draws y,Sigma,mu,tau,sigma,Lambda,e
#                                R. McCulloch 12/28/04
# 
# arguments:
#    Data:
#     all components are required:
#       k:  integer giving the scale of the responses, each observation is an integer from 1,2,...k
#       x:  data, num rows=number of respondents, num columns = number of questions
#    Prior:
#     all components are optional
#       nu,V: Sigma ~ IW(nu,V)
#       mubar,Am: mu ~N(mubar,Am^{-1})
#       gsigma: grid for sigma
#       gl11,gl22,gl12: grids for ij element of Lamda
#       Lambdanu,LambdaV: Lambda ~ IW(Lambdanu,LambdaV)
#       ge: grid for e
#    Mcmc:
#     all components are optional (but you would typically want to specify R= number of draws)
#       R: number of mcmc iterations
#       keep: frequency with which draw is kept
#       ndghk: number of draws for ghk
#       printevery: how often to print out how many draws are done
#       e,y,mu,Sigma,sigma,tau,Lamda: initial values for the state
#       doe, ...doLambda: indicates whether draw should be made
# output:
#    List with draws of each of Sigma,mu,tau,sigma,Lambda,e
#    eg. result$Sigma is the draws of Sigma
#    Each component is a matrix expept e, which is a vector
#    for the matrices Sigma and Lambda each row transpose of the Vec
#    eg. result$Lambda has rows (Lambda11,Lamda21,Lamda12,Lamda22)

#
# define functions needed
#
# -----------------------------------------------------------------------------------
rlpx = function(x,e,k,mu,tau,Sigma,sigma,nd=500) {
n=nrow(x); p = ncol(x)
cc = cgetC(e,k)
L=t(chol(Sigma))
lpv = rep(0,n)
offset = p*log(k)
for(i in 1:n) {
   Li = sigma[i]*L
   a = cc[x[i,]]-mu-tau[i]; b = cc[x[i,]+1]-mu-tau[i]
   ghkres = rghk(Li,a,b,nd)
   lghkres = log(ghkres)
   if(is.nan(lghkres)) {
      #print("nan in ghk:")
      #print(paste('ghkres: ',ghkres))
      lghkres = log(1e-320)
   }
   if(is.infinite(lghkres)) {
      #print("infinite in ghk:")
      #print(paste('ghkres: ',ghkres))
      lghkres = log(1e-320)
   }
   lpv[i] = lghkres + offset
}
sum(lpv)
}
rghk = function(L,a,b,nd) {
.C('ghk',as.double(L),as.double(a),as.double(b),as.integer(nrow(L)),
          as.integer(nd),res=double(1))$res
}
condd = function(Sigma) {
p = nrow(Sigma)
Si = solve(Sigma)
cbeta = matrix(0,p-1,p)
for(i in 1:p) {
ind = (1:p)[-i]
cbeta[,i] = -Si[ind,i]/Si[i,i]
}
list(beta=cbeta,s=sqrt(1/diag(Si)))
}
pandterm = function(message) { stop(paste("in rscaleUsage: ",message),call.=FALSE) }
myin = function(i,ind) {i %in% ind}
getS = function(Lam,n,moms) {
S=matrix(0.0,2,2)
S[1,1] = (n-1)*moms[3] + n*moms[1]^2
S[1,2] = (n-1)*moms[4] + n*moms[1]*(moms[2]-Lam[2,2])
S[2,1] = S[1,2]
S[2,2] = (n-1)*moms[5] + n*(moms[2]-Lam[2,2])^2
S
}
llL = function(Lam,n,S,V,nu) {
dlam = Lam[1,1]*Lam[2,2]-Lam[1,2]^2
M = (S+V) %*%  chol2inv(chol(Lam))
ll = -.5*(n+nu+3)*log(dlam) -.5*sum(diag(M))
}
ispd = function(mat,d=nrow(mat)) {
if(!is.matrix(mat)) {
res = FALSE
} else if(!((nrow(mat)==d) & (ncol(mat)==d))) {
res = FALSE
} else {
diff = (t(mat)+mat)/2 - mat
perdiff = sum(diff^2)/sum(mat^2)
res = ((det(mat)>0) & (perdiff < 1e-10))
}
res
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# print out components of inputs ----------------------------------------------
cat('\nIn function rscaleUsage\n\n')
if(!missing(Data)) {
cat('   Data has components: ')
cat(paste(names(Data),collapse=' ')[1],'\n')
}
if(!missing(Prior)) {
cat('   Prior has components: ')
cat(paste(names(Prior),collapse=' ')[1],'\n')
}
if(!missing(Mcmc)) {
cat('   Mcmc has components: ')
cat(paste(names(Mcmc),collapse=' ')[1],'\n')
}
cat('\n')
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# process Data argument --------------------------
if(missing(Data)) {pandterm("Requires Data argument - list of k=question scale and x = data")}
if(is.null(Data$k)) {
   pandterm("k not specified")
} else {
   k = as.integer(Data$k)
   if(!((k>0) & (k<50))) {pandterm("Data$k must be integer between 1 and 50")}
}
if(is.null(Data$x)) {
   pandterm('x (the data), not specified')
} else {
   if(!is.matrix(Data$x)) {pandterm('Data$x must be a matrix')}
   x = matrix(as.integer(Data$x),nrow=nrow(Data$x))
   checkx = sum(sapply(as.vector(x),myin,1:k))
   if(!(checkx == nrow(x)*ncol(x))) {pandterm('each element of Data$x must be in 1,2...k')}
   p = ncol(x)
   n = nrow(x)
   if((p<2) | (n<1)) {pandterm(paste('invalid dimensions for x: nrow,ncol: ',n,p))}
}
# ++++++++++++++++++++++++++++++++++++++++++++++++

# process Mcmc argument ---------------------

#run mcmc
R = as.integer(1000)
keep = as.integer(1)
ndghk= as.integer(100)
printevery = as.integer(10)
if(!missing(Mcmc)) {
if(!is.null(Mcmc$R))              { R = as.integer(Mcmc$R) }
if(!is.null(Mcmc$keep))           { keep = as.integer(Mcmc$keep) }
if(!is.null(Mcmc$ndghk))          { ndghk = as.integer(Mcmc$ndghk) }
if(!is.null(Mcmc$printevery))     { printevery = as.integer(Mcmc$printevery) }
}
if(R<1) { pandterm('R must be positive')}
if(keep<1) { pandterm('keep must be positive') }
if(ndghk<1) { pandterm('ndghk must be positive') }
if(printevery<1) { pandterm('printevery must be positive') }


#state
y = matrix(as.double(x),nrow=nrow(x))
mu = apply(y,2,mean)
Sigma = var(y)
tau = rep(0,n)
sigma = rep(1,n)
#Lamda = matrix(c(3.7,-.22,-.22,.32),ncol=2)
#Lamda = matrix(c((k/4)^2,(k/4)*.5*(-.2),0,.25),nrow=2); Lamda[1,2]=Lamda[2,1]
Lamda = matrix(c(4,0,0,.5),ncol=2)
e=0
if(!missing(Mcmc)) {
if(!is.null(Mcmc$y))         { y = Mcmc$y }
if(!is.null(Mcmc$mu))        { mu = Mcmc$mu }
if(!is.null(Mcmc$Sigma))     { Sigma = Mcmc$Sigma }
if(!is.null(Mcmc$tau))       { tau = Mcmc$tau }
if(!is.null(Mcmc$sigma))     { sigma = Mcmc$sigma }
if(!is.null(Mcmc$Lambda))    { Lamda = Mcmc$Lambda }
if(!is.null(Mcmc$e))         { e = Mcmc$e }
}
if(!ispd(Sigma,p)) { pandterm(paste('Sigma must be positive definite with dimension ',p)) }
if(!ispd(Lamda,2)) { pandterm(paste('Lambda must be positive definite with dimension ',2)) }
if(!is.vector(mu)) { pandterm('mu must be a vector') }
if(length(mu) != p) { pandterm(paste('mu must have length ',p)) }
if(!is.vector(tau)) { pandterm('tau must be a vector') }
if(length(tau) != n) { pandterm(paste('tau must have length ',n)) }
if(!is.vector(sigma)) { pandterm('sigma must be a vector') }
if(length(sigma) != n) { pandterm(paste('sigma must have length ',n)) }
if(!is.matrix(y)) { pandterm('y must be a matrix') }
if(nrow(y) != n) { pandterm(paste('y must have',n,'rows')) }
if(ncol(y) != p) { pandterm(paste('y must have',p,'columns')) }

#do draws
domu=TRUE
doSigma=TRUE
dosigma=TRUE
dotau=TRUE
doLamda=TRUE
doe=TRUE
if(!missing(Mcmc)) {
if(!is.null(Mcmc$domu))        { domu = Mcmc$domu }
if(!is.null(Mcmc$doSigma))     { doSigma = Mcmc$doSigma }
if(!is.null(Mcmc$dotau))       { dotau = Mcmc$dotau }
if(!is.null(Mcmc$dosigma))     { dosigma = Mcmc$dosigma }
if(!is.null(Mcmc$doLambda))    { doLamda = Mcmc$doLambda }
if(!is.null(Mcmc$doe))         { doe = Mcmc$doe }
}


#++++++++++++++++++++++++++++++++++++++

#process Prior argument ----------------------------------
nu = p+3
V= nu*diag(p)
mubar = matrix(rep(k/2,p),ncol=1)
Am = .0001*diag(p)
gs = 500
gsigma = 6*(1:gs)/gs
gl11 = .1 + 5.9*(1:gs)/gs
gl22 = .1 + 2.0*(1:gs)/gs
#gl12 = -.8 + 1.6*(1:gs)/gs
gl12 = -2.0 + 4*(1:gs)/gs
nuL=20
VL = (nuL-3)*Lamda
ge = -.1+.2*(0:gs)/gs

if(!missing(Prior)) {
if(!is.null(Prior$nu))       { nu = Prior$nu; V = nu*diag(p) }
if(!is.null(Prior$V))        { V = Prior$V }
if(!is.null(Prior$mubar))    { mubar = matrix(Prior$mubar,ncol=1) }
if(!is.null(Prior$Am))       { Am = Prior$Am }
if(!is.null(Prior$gsigma))   { gsigma = Prior$gsigma }
if(!is.null(Prior$gl11))     { gl11 = Prior$gl11 }
if(!is.null(Prior$gl22))     { gl22 = Prior$gl22 }
if(!is.null(Prior$gl12))     { gl12 = Prior$gl12 }
if(!is.null(Prior$Lambdanu)) { nuL = Prior$Lambdanu; VL = (nuL-3)*Lamda }
if(!is.null(Prior$LambdaV))  { VL = Prior$LambdaV }
if(!is.null(Prior$ge))       { ge = Prior$ge }
}
if(!ispd(V,p)) { pandterm(paste('V must be positive definite with dimension ',p)) }
if(!ispd(Am,p)) { pandterm(paste('Am must be positive definite with dimension ',p)) }
if(!ispd(VL,2)) { pandterm(paste('VL must be positive definite with dimension ',2)) }
if(nrow(mubar) != p) { pandterm(paste('mubar must have length',p)) }
#++++++++++++++++++++++++++++++++++++++++

#print out run info -------------------------
cat('   n,p,k: ', n,p,k,'\n')
cat('   R,keep,ndghk,printevery: ', R,keep,ndghk,printevery,'\n')
cat('\n')
cat('   Data:\n')
cat('      x11,n1,1p,np: ',x[1,1],x[n,1],x[1,p],x[n,p],'\n\n')
cat('   Prior:\n')
cat('      ','nu: ',nu,'\n')
cat('      ','V11,pp/nu: ',V[1,1]/nu,V[p,p]/nu,'\n')
cat('      ','mubar1,p: ',mubar[1],mubar[p],'\n')
cat('      ','Am11,pp: ',Am[1,1],Am[p,p],'\n')
cat('      ','Lambdanu: ',nuL,'\n')
cat('      ','LambdaV11,22/(Lambdanu-3): ',VL[1,1]/(nuL-3),VL[2,2]/(nuL-3),'\n')
cat('      ','mubar1,p: ',mubar[1],mubar[p],'\n')
cat('      ','sigma grid, 1,',length(gsigma),': ',gsigma[1],', ',gsigma[length(gsigma)],'\n')
cat('      ','Lambda11 grid, 1,',length(gl11),': ',gl11[1],', ',gl11[length(gl11)],'\n')
cat('      ','Lambda12 grid, 1,',length(gl12),': ',gl12[1],', ',gl12[length(gl12)],'\n')
cat('      ','Lambda22 grid, 1,',length(gl22),': ',gl22[1],', ',gl22[length(gl22)],'\n')
cat('      ','e grid, 1,',length(ge),': ',ge[1],', ',ge[length(ge)],'\n')
cat('      ','draw e: ',doe,'\n')
cat('      ','draw Lambda: ',doLamda,'\n')
#++++++++++++++++++++++++++++++++++++++++++++

nk = floor(R/keep)
ndpost = nk*keep
drSigma=matrix(0.0,nk,p^2)
drmu = matrix(0.0,nk,p)
drtau = matrix(0.0,nk,n)
drsigma = matrix(0.0,nk,n)
drLamda = matrix(0.0,nk,4)
dre = rep(0,nk)

itime = proc.time()[3]
cat("Mcmc Iteration (est time to end - min)",'\n')
for(rep in 1:ndpost) {
   if(1) { # y
      cc = cgetC(e,k)
      bs = condd(Sigma)
      y = matrix(.C('dy',as.integer(p),as.integer(n),y=as.double(t(y)),as.integer(t(x)),as.double(cc),as.double(mu),as.double(bs$beta),as.double(bs$s),
	                                        as.double(tau),as.double(sigma))$y,ncol=p,byrow=TRUE)
   }
   if(doSigma) { #Sigma
      Res = (t(t(y)-mu)-tau)/sigma
      S = crossprod(Res)
      Sigma = rwishart(nu+n,chol2inv(chol(V+S)))$IW
   }
   if(domu) { #mu
      yd = y-tau
      Si = chol2inv(chol(Sigma))
      Vmi = sum(1/sigma^2)*Si + Am
      R = chol(Vmi)
      Ri = backsolve(R,diag(rep(1,p)))
      Vm = chol2inv(chol(Vmi))
      mm = Vm %*% (Si %*% (t(yd) %*% matrix(1/sigma^2,ncol=1)) + Am %*% mubar)
      mu = as.vector(mm + Ri %*% matrix(rnorm(p),ncol=1))
   }
   if(dotau) { #tau
      Ai = Lamda[1,1] - (Lamda[1,2]^2)/Lamda[2,2]
      A = 1.0/Ai
      onev = matrix(1.0,p,1)
      R = chol(Sigma)
      xx = backsolve(R,onev,transpose=TRUE)
      yy = backsolve(R,t(y)-mu,transpose=TRUE)
      xtx = sum(xx^2)
      xty = as.vector(t(xx) %*% yy)
      beta = A*Lamda[1,2]/Lamda[2,2]
      for(j in 1:n) {
	 s2 = xtx/sigma[j]^2   + A
         s2 = 1.0/s2 
	 m = s2*((xty[j]/sigma[j]^2) + beta*(log(sigma[j])-Lamda[2,2]))
	 tau[j] = m + sqrt(s2)*rnorm(1)
      }
   }
   if(dosigma) { #sigma
      R = chol(Sigma)
      eps = backsolve(R,t(y-tau)-mu,transpose=TRUE)
      ete = as.vector(matrix(rep(1,p),nrow=1) %*% eps^2)
      a= Lamda[2,2]
      b= Lamda[1,2]/Lamda[1,1]
      s=sqrt(Lamda[2,2]-(Lamda[1,2]^2/Lamda[1,1]))
      for(j in 1:n) {
	 pv = -(p+1)*log(gsigma) -.5*ete[j]/gsigma^2 -.5*((log(gsigma)-(a+b*tau[j]))/s)^2
	 pv = exp(pv-max(pv))
	 pv = pv/sum(pv)
	 sigma[j] = sample(gsigma,size=1,prob=pv)
      }
   }
   
   if(doLamda) { # Lamda
      h=log(sigma)
      dat = cbind(tau,h)
      temp = var(dat)
      moms = c(mean(tau),mean(h),temp[1,1],temp[1,2],temp[2,2])

      SS = getS(Lamda,n,moms)
      rgl11 = gl11[gl11 > (Lamda[1,2]^2/Lamda[2,2])]
      ng = length(rgl11)
      pv = rep(0,ng)
      for(j in 1:ng) {
         Lamda[1,1] = rgl11[j]
         pv[j] = llL(Lamda,n,SS,VL,nuL)
      }
      pv = exp(pv-max(pv)); pv = pv/sum(pv)
      Lamda[1,1] = sample(rgl11,size=1,prob=pv)

      rgl12 = gl12[(gl12<sqrt(Lamda[1,1]*Lamda[2,2])) & (gl12>-sqrt(Lamda[1,1]*Lamda[2,2]))]
      ng = length(rgl12)
      pv = rep(0,ng)
      for(j in 1:ng) {
         Lamda[1,2] = rgl12[j]; Lamda[2,1]=Lamda[1,2]
         pv[j] = llL(Lamda,n,SS,VL,nuL)
      }
      pv = exp(pv-max(pv)); pv = pv/sum(pv)
      Lamda[1,2] = sample(rgl12,size=1,prob=pv)
      Lamda[2,1]=Lamda[1,2]

      rgl22 = gl22[gl22 > (Lamda[1,2]^2/Lamda[1,1])]
      ng = length(rgl22)
      pv = rep(0,ng)
      for(j in 1:ng) {
         Lamda[2,2] = rgl22[j]
	 SS = getS(Lamda,n,moms)
         pv[j] = llL(Lamda,n,SS,VL,nuL)
      }
      pv = exp(pv-max(pv)); pv = pv/sum(pv)
      Lamda[2,2] = sample(rgl22,size=1,prob=pv)
   }

   if(doe) { # e
      ng = length(ge)
      ei = which.min(abs(e-ge))
      if(ei==1) {
         pi =2
         qr = .5
      } else if(ei==ng) {
         pi = ng-1
         qr = .5
      } else {
         pi = ei + rbinom(1,1,.5)*2-1
         qr = 1
      }
      eold = ge[ei]
      eprop = ge[pi]
      llold = rlpx(x,eold,k,mu,tau,Sigma,sigma,ndghk)
      llprop = rlpx(x,eprop,k,mu,tau,Sigma,sigma,ndghk)
      lrat = llprop - llold + log(qr)
      if(lrat>0) {
         e = eprop
      } else {
         paccept = min(1,exp(lrat))
         e = ifelse(rbinom(1,1,paccept),eprop,eold)
      }
   }
   mkeep = rep/keep
   if(mkeep == floor(mkeep)) {
      drSigma[mkeep,] = Sigma
      drmu[mkeep,] = mu
      drtau[mkeep,] = tau
      drsigma[mkeep,] = sigma 
      drLamda[mkeep,] = Lamda
      dre[mkeep] = e
   }
   if((rep/printevery)==floor(rep/printevery)) {
      ctime = proc.time()[3]
      timetoend = ((ctime-itime)/rep)*(ndpost-rep)
      cat(rep,' (', round(timetoend/60,1), ') \n')
      fsh()
   }
}
ctime = proc.time()[3]
cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')
return(list(Sigmadraw=drSigma,mudraw=drmu,taudraw = drtau,
 sigmadraw=drsigma,Lambdadraw=drLamda,edraw=dre))
}




