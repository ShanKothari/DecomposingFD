## calculating structured gamma diversity
## currently non-functional
FTD.gamma.str<-function(tdmat,spmat,abund=F,q=1){
  St<-sum(spmat>0)
  nsp.comm<-rowSums(spmat>0)
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
  
  sp.rep<-which(spmat>0,arr.ind=T)[,2]
  td.sp<-tdmat[sp.rep,sp.rep]
  
  ## is this the correct formula for weights?
  ## or should proportional abundances in less speciose communities count for more?
  ## alternative:
  ## weights<-spmat[which(spmat>0)]/n.comm
  spmat.sp<-spmat*nsp.comm
  weights<-spmat.sp[which(spmat.sp>0)]/St
  
  td.norm<-diag(weights) %*% td.sp %*% diag(weights)
  
  M.gamma<-sum(td.norm)
  M.gamma.prime<-M.gamma*St/(St-1)

  fij<-td.norm/M.gamma
  if(q==1){
    fijlog<-fij*log(fij)
    fijlog[is.na(fijlog)]<-0
    Ht.gamma<-exp(-1*sum(fijlog))
  } else {
    Ht.gamma<-sum(fij^q)^(1/(1-q))
  }
  
  qDT.gamma<-(1+sqrt(1+4*Ht.gamma))/2
  qDTM.gamma<-1+qDT.gamma*M.gamma
  Et.gamma<-qDT.gamma/St
  
  list(St=St,M.gamma=M.gamma,Ht.gamma=Ht.gamma,Et.gamma=Et.gamma,qDT.gamma=qDT.gamma,qDTM.gamma=qDTM.gamma)
}