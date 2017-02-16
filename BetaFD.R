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

comm.disp.mat<-function(tdmat,spmat,weighted=FALSE){
  n.comm<-nrow(spmat)
  disp.mat<-outer(1:n.comm,1:n.comm,FUN=Vectorize(function(i,j) comm.disp(tdmat,com1=spmat[i,],com2=spmat[j,])))
  if(weighted==T){
    nsp.comm<-rowSums(spmat>0)
    disp.mat.weight<-diag(nsp.comm) %*% disp.mat %*% diag(nsp.comm)
    return(disp.mat.weight)
  } else {
    return(disp.mat)
  }
}

M.gamma.pw<-function(tdmat,spmat){
  M.gamma<-function(tdmat,com1,com2){
    c.ind<-c(which(com1>0),which(com2>0))
    M.c<-mean(tdmat[c.ind,c.ind])/length(c.ind)^2
    return(M.c)
  }
  n.comm<-nrow(spmat)
  M.gamma.mat<-outer(1:n.comm,1:n.comm,function(i,j) M.gamma(tdmat,spmat[i,],spmat[j,]))
  return(M.gamma.mat)
}

M.beta.pairwise<-function(tdmat,spmat,norm=F){
  nsp.comm<-rowSums(spmat>0)
  n.comm<-nrow(spmat)
  nsp.pair<-outer(1:n.comm,1:n.comm,function(i,j) nsp.comm[i]+nsp.comm[j])
  disp.mat.weight<-comm.disp.mat(tdmat,spmat,weighted=T)
  M.beta.pw<-2*disp.mat.weight/nsp.pair^2
  if(norm==T){
    M.beta.norm<-M.beta.pw/M.gamma.pw(tdmat,spmat)
    return(M.beta.norm)
  } else {
    return(M.beta.pw)
  }
}

FTD.beta<-function(tdmat,spmat,q=1){
  n.comm<-nrow(spmat)
  disp.mat.weight<-comm.disp.mat(tdmat,spmat,weighted=T)
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