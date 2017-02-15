## measuring the dispersion between two communities
## based on presence/absence
comm.disp<-function(tdmat,com1,com2){
  ## mean distance between species in community 1 and 2
  mAB<-mean(tdmat[com1>0,com2>0])
  ## correction so that dispersion=0 if communities have the same species
  mAA<-mean(tdmat[com1>0,com1>0])
  mBB<-mean(tdmat[com2>0,com2>0])
  dmAB<-mAB-(mAA+mBB)/2
  return(dmAB)
}

FTD.beta<-function(tdmat,spmat,q=1){
  n.comm<-nrow(spmat)
  disp.mat<-outer(1:n.comm,1:n.comm,FUN=Vectorize(function(i,j) comm.disp(tdmat,com1=spmat[i,],com2=spmat[j,])))

  nsp.comm<-rowSums(spmat>0)
  disp.mat.weight<-diag(nsp.comm) %*% disp.mat %*% diag(nsp.comm)
  M.beta<-sum(disp.mat.weight)/sum(spmat>0)^2
  
  fAB<-disp.mat/sum(disp.mat)
  if(q==1){
    fABlog<-fAB*log(fAB)
    fABlog[is.na(fABlog)]<-0
    Ht.beta<-exp(-1*sum(fABlog))
  } else {
    Ht.beta<-sum(fAB^q)^(1/(1-q))
  }
  qDT.beta<-(1+sqrt(1+4*Ht.beta))/2
  qDTM.beta<-1+qDT.beta*M.beta
  Et.beta<-qDT.beta/n.comm
  
  list(n.comm=n.comm,q=q,M.beta=M.beta,Ht.beta=Ht.beta,qDT.beta=qDT.beta,qDTM.beta=qDTM.beta,disp.mat=disp.mat)
}