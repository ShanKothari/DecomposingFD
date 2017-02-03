## in progress!
## beginnings of scripts to calculate Scheiner's functional beta-diversity

## measuring the dispersion between two communities
## based on presence/absence
comm.disp<-function(tdmat,com1,com2){
  ## mean distance between species in community 1 and 2
  mAB<-mean(tdmat[com1>0,com2>0])
  ## correction so that dispersion=0 if communities have the same species
  mAA<-mean(tdmat[com1>0,com1>0])
  mBB<-mean(tdmat[com2>0,com2>0])
  nsp.com1<-sum(com1>0)
  nsp.com2<-sum(com2>0)
  dmAB<-mAB-(nsp.com1*mAA+nsp.com2*mBB)/(nsp.com1+nsp.com2)
  dmAB.weight<-dmAB*nsp.com1*nsp.com2
  disp<-list(dmAB,dmAB.weight)
  return(disp)
}

## note: this occasionally yields negative values
## even when using distances like Euclidean or Manhattan distance
## Scheiner et al. claim this is impossible (Appendix, Dispersion of Communities)
## due to concavity -- I suspect that it is, in fact, possible
## and concavity does not apply

## measuring the dispersion among a set of communities
## wraps the above across a communities x species matrix
multicomm.disp<-function(tdmat,spmat){
  disp.mat<-outer(1:nrow(spmat),1:nrow(spmat),FUN=Vectorize(function(i,j) comm.disp(tdmat,com1=spmat[i,],com2=spmat[j,])[[2]]))
  MB<-(sum(disp.mat)/sum(spmat>0)^2)*nrow(spmat)/(nrow(spmat)-1)
  disp.list<-list(disp.mat,MB)
  return(disp.list)
}