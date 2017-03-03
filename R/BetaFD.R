## measuring the dispersion between two communities
## based on presence/absence
comm.disp<-function(tdmat,com1,com2){
  
  if(identical(sum(com1),1)==F || identical(sum(com2),1)==F){
    com1<-com1/sum(com1)
    com2<-com2/sum(com2)
    warning("input proportional abundances do not sum to 1; summation to 1 forced")
  }
  
  ## mean distance between species in community 1 and 2
  mAB<-sum(diag(com1) %*% tdmat %*% diag(com2))
  ## correction so that dispersion=0 if communities have the same species
  mAA<-sum(diag(com1) %*% tdmat %*% diag(com1))
  mBB<-sum(diag(com2) %*% tdmat %*% diag(com2))
  dmAB<-mAB-(mAA+mBB)/2
  return(dmAB)
}


comm.disp.mat<-function(tdmat,spmat,abund=F,sp.weighted=FALSE){
  n.comm<-nrow(spmat)
  
  if(abund==F){
    spmat[spmat>0]<- 1
    spmat<-spmat/rowSums(spmat)
  }
  
  ## if any row doesn't sum to 1, coerce summation
  if(FALSE %in% sapply(rowSums(spmat),function(x) identical(x,1))){
    spmat<-spmat/rowSums(spmat)
    warning("proportional abundances don't always sum to 1; summation to 1 forced")
  }
  
  disp.mat<-outer(1:n.comm,1:n.comm,FUN=Vectorize(function(i,j) comm.disp(tdmat,com1=spmat[i,],com2=spmat[j,])))
  if(sp.weighted==T){
    nsp.comm<-rowSums(spmat>0)
    disp.mat.weight<-diag(nsp.comm) %*% disp.mat %*% diag(nsp.comm)
    return(disp.mat.weight)
  } else {
    return(disp.mat)
  }
}

M.gamma.pairwise<-function(tdmat,spmat,abund=F){
  
  if(abund==F){
    spmat[spmat>0]<- 1
    spmat<-spmat/rowSums(spmat)
  }
  
  ## if any row doesn't sum to 1, coerce summation
  if(FALSE %in% sapply(rowSums(spmat),function(x) identical(x,1))){
    spmat<-spmat/rowSums(spmat)
    warning("proportional abundances don't always sum to 1; summation to 1 forced")
  }
  
  M.gamma<-function(tdmat,com1,com2){
    nsp1<-sum(com1>0)
    nsp2<-sum(com2>0)
    c.ind<-c(which(com1>0),which(com2>0))
    ## is this the correct abundance-weighted M.gamma?
    ## or should proportional abundances in less speciose communities count for more?
    ## alternative:
    ## c.abund<-c(com1[which(com1>0)],com2[which(com2>0)])/2
    c.abund<-c(com1[which(com1>0)]*nsp1,com2[which(com2>0)]*nsp2)/(nsp1+nsp2)
    M.c<-sum(c.abund %*% tdmat[c.ind,c.ind] %*% c.abund)
    return(M.c)
  }
  
  n.comm<-nrow(spmat)
  M.gamma.mat<-outer(1:n.comm,1:n.comm,FUN=Vectorize(function(i,j) M.gamma(tdmat,spmat[i,],spmat[j,])))
  return(M.gamma.mat)
}

M.beta.pairwise<-function(tdmat,spmat,abund=F,norm=F){
  nsp.comm<-rowSums(spmat>0)
  n.comm<-nrow(spmat)
  nsp.pair<-outer(1:n.comm,1:n.comm,function(i,j) nsp.comm[i]+nsp.comm[j])
  disp.mat.weight<-comm.disp.mat(tdmat,spmat,abund=abund,sp.weighted=T)
  ## factor of 2 to include distance from A to B and B to A
  M.beta.mat<-2*disp.mat.weight/nsp.pair^2
  if(norm==T){
    M.beta.norm<-M.beta.mat/M.gamma.pairwise(tdmat,abund=abund,spmat)
    M.beta.norm[is.na(M.beta.norm)]<-0
    return(M.beta.norm)
  } else {
    return(M.beta.mat)
  }
}

FTD.beta<-function(tdmat,spmat,abund=F,q=1){
  nsp<-sum(colSums(spmat)>0)
  St<-sum(spmat>0)
  n.comm<-nrow(spmat)
  
  disp.mat.weight<-comm.disp.mat(tdmat,spmat,abund=abund,sp.weighted=T)
  M.beta<-sum(disp.mat.weight)/St^2
  M.beta.prime<-M.beta*n.comm/(n.comm-1)
  
  fAB<-disp.mat.weight/sum(disp.mat.weight)
  if(q==1){
    fABlog<-fAB*log(fAB)
    fABlog[is.na(fABlog)]<-0
    Ht.beta<-exp(-1*sum(fABlog))
  } else if(q==0){
    Ht.beta<-sum(fAB>0)
  } else {
    Ht.beta<-sum(fAB^q)^(1/(1-q))
  }
  qDT.beta<-(1+sqrt(1+4*Ht.beta))/2
  qDTM.beta<-1+qDT.beta*M.beta
  Et.beta<-qDT.beta/n.comm
  
  list(nsp=nsp,St=St,n.comm=n.comm,q=q,M.beta=M.beta,M.beta.prime=M.beta.prime,Ht.beta=Ht.beta,qDT.beta=qDT.beta,qDTM.beta=qDTM.beta,disp.mat.weight=disp.mat.weight)
}