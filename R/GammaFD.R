## method 1
FTD.gamma.str<-function(tdmat,spmat,q=1){
  St<-sum(spmat>0)
  sp.rep<-which(spmat>0,arr.ind=T)[,2]
  td.sp<-tdmat[sp.rep,sp.rep]
  weights<-rep(1/St,St)
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