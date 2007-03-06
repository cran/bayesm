plot.bayesm.nmix=function(x,names,burnin=trunc(.1*nrow(probdraw)),Grid,bi.sel,nstd=2,...){
#
# S3 method to plot normal mixture marginal and bivariate densities
#     nmixlist is a list of 3 components, nmixlist[[1]]: array of mix comp prob draws,
#     mmixlist[[2]] is not used, nmixlist[[3]] is list of draws of components
#     P. Rossi 2/08
#
  nmixlist=x
  if(mode(nmixlist) != "list") stop(" Argument must be a list \n")
  probdraw=nmixlist[[1]]; compdraw=nmixlist[[3]]
  if(!is.matrix(probdraw)) stop(" First element of list (probdraw) must be a matrix \n")
  if(mode(compdraw) != "list") stop(" Third element of list (compdraw) must be a list \n")
  op=par(no.readonly=TRUE)
  on.exit(par(op))
  R=nrow(probdraw)
  if(R < 100) {cat(" fewer than 100 draws submitted \n"); return(invisible())}
  datad=length(compdraw[[1]][[1]]$mu)
  if(missing(bi.sel)) bi.sel=list(c(1,2))  # default to the first pair of variables
  ind=as.integer(seq(from=(burnin+1),to=R,length.out=max(200,trunc(.05*R))))
  if(missing(Grid)){
     out=momMix(probdraw[ind,],compdraw[ind])
     mu=out$mu
     sd=out$sd
     Grid=matrix(0,nrow=50,ncol=datad)
     for(i in 1:datad ) Grid[,i]=c(seq(from=(mu[i]-nstd*sd[i]),to=(mu[i]+nstd*sd[i]),length=50))
     }
  #
  #  plot posterior mean of marginal densities
  #
  mden=eMixMargDen(Grid,probdraw[ind,],compdraw[ind])
  nx=datad
  if(nx==1) par(mfrow=c(1,1)) 
  if(nx==2) par(mfrow=c(2,1))
  if(nx==3) par(mfrow=c(3,1))
  if(nx==4) par(mfrow=c(2,2))
  if(nx>=5) par(mfrow=c(3,2))

  if(missing(names)) {names=as.character(1:nx)}
  for(index in 1:nx){
        if(index == 2) par(ask=dev.interactive)
        plot(range(Grid[,index]),c(0,1.1*max(mden[,index])),type="n",xlab="",ylab="density")
        title(names[index])
        lines(Grid[,index],mden[,index],col="black",lwd=1.5)
        polygon(c(Grid[1,index],Grid[,index],Grid[nrow(Grid),index]),c(0,mden[,index],0),col="magenta")
  }
  #
  # now plot bivariates in list bi.sel
  #
  par(ask=dev.interactive())
  nsel=length(bi.sel)
  ngrid=50
  den=array(0,dim=c(nsel,ngrid,ngrid))
  lstxixj=NULL
  for(sel in 1:nsel){
      i=bi.sel[[sel]][1]
      j=bi.sel[[sel]][2]
      rxi=range(Grid[,i])
      rxj=range(Grid[,j])
      xi=seq(from=rxi[1],to=rxi[2],length.out=ngrid)
      xj=seq(from=rxj[1],to=rxj[2],length.out=ngrid)
      lstxixj[[sel]]=list(xi,xj)
      for(elt in ind){
         den[sel,,]=den[sel,,]+mixDenBi(i,j,xi,xj,probdraw[elt,],compdraw[[elt]])
      }
   }     
  nx=nsel
  par(mfrow=c(1,1))

  for(index in 1:nx){
        xi=unlist(lstxixj[[index]][1])
        xj=unlist(lstxixj[[index]][2])
        xlabtxt=paste("Var ",bi.sel[[index]][1],sep="")
        ylabtxt=paste("Var ",bi.sel[[index]][2],sep="")
        image(xi,xj,den[index,,],col=terrain.colors(100),xlab=xlabtxt,ylab=ylabtxt)
        contour(xi,xj,den[index,,],add=TRUE,drawlabels=FALSE)
  }

  invisible()
}

