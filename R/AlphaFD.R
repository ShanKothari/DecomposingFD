## FTD: A function to calculate Scheiner et al. (2016)'s functional
## trait dispersion, written by S.A. Kothari
## with a few lines of code shamelessly stolen from J.J. Grossman

## contains two functions:
## 1. FTD uses a trait distance matrix to calculate functional diversity
## 2. FTD.comm applies FTD across the rows of a community x species matrix
## to calculate FTD for each community

## FTD requires:
## matrix or dist containing species functional distances (tdmat)
## tdmat should be rescaled to 0<=dij<=1
## Hill coefficienct q
## weights are proportions of each species in each plot

FTD<-function(tdmat,q=1,abund=F,weights=NULL){
  ## contingency for one-species communities
  if(length(tdmat)==1 && tdmat==0){
    tdmat<-as.matrix(tdmat)
  }
  ## is the input a matrix or dist? if not...
  if(!(class(tdmat) %in% c("matrix","dist"))){
    stop("distances must be class dist or class matrix")
  } else if(class(tdmat)=="matrix" && !isSymmetric(tdmat)){
    warning("matrix not symmetric")
  } else if(class(tdmat)=="dist"){
    tdmat<-as.matrix(tdmat)
  }
  
  if(max(tdmat)>1 || min(tdmat)<0){
    tdmat<-(tdmat-min(tdmat))/(max(tdmat)-min(tdmat))
    warning("trait distances must be between 0 and 1; rescaling")
  }
  
  nsp<-nrow(tdmat)
  ## if abund=T but no weights are provided, abundances are assumed equal
  if(is.null(weights)){
    weights<-rep(1/nsp,nsp)
    if(abund==T){
      warning("no weights provided; weights assumed equal")
    }
  }
  
  if(identical(sum(weights),1)==F){
    weights<-weights/sum(weights)
    warning("input proportional abundances do not sum to 1; summation to 1 forced")
  }
  
  ## if abund=F, tdmat.abund=dij
  ## if abund=T, tdmat.abund is weighted by proportional abundance of i and j
  tdmat.abund<-diag(weights^abund) %*% tdmat %*% diag(weights^abund)
  dijsum<-sum(tdmat.abund)
  ## calculates Scheiner's M (if abund=F) or Rao's Q (if abund=T)
  M<-dijsum/nsp^(2*(1-abund))
  M.prime<-M*nsp/(nsp-1)
  fij<-tdmat.abund/dijsum
  
  ## calculating qH
  ## fork -- if q=1, 1/(1-q) is undefined, so we use an analogue
  ## of the formula for Shannon-Weiner diversity
  ## if q!=1, we can calculate explicitly
  if(q==1){
    fijlog<-fij*log(fij)
    fijlog[is.na(fijlog)]<-0
    Ht<-exp(-1*sum(fijlog))
  } else {
    Ht<-sum(fij^q)^(1/(1-q))
  }
  
  ## getting qDT, qDTM, and Et from Ht
  qDT<-(1+sqrt(1+4*Ht))/2
  qDTM<-1+qDT*M
  Et<-qDT/nsp
  
  list(nsp=nsp,q=q,M=M,M.prime=M.prime,Ht=Ht,Et=Et,qDT=qDT,qDTM=qDTM)
}

## wrapper for the above function across multiple communities
## requires distance matrix w/ all species, scaled to 0-1
## community data matrix (communities x species)
## with species in same order as distance matrix
## if abund=T, values in community data matrix are treated as abundances
## if match.names=T, the code will match species names across the
## trait distance matrix and comm data matrix and rearrange the latter

FTD.comm<-function(tdmat,spmat,q=1,abund=F,match.names=F){
  
  n.comm<-nrow(spmat)
  if(match.names==T){
    sp.arr<-match(rownames(as.matrix(tdmat)),colnames(spmat))
    spmat<-spmat[,sp.arr]
  }
  
  ## define an internal function to apply FTD to a subsetted distance matrix
  select.FTD<-function(tdmat,spvec,q,abund){
    sp<-which(as.logical(spvec))
    tdmat.comm<-as.matrix(tdmat)[sp,sp]
    ## automatically converts abundances to proportions
    output<-FTD(tdmat=tdmat.comm,q=q,abund=abund,weights=as.numeric(spvec[sp])/sum(spvec[sp]))
    output<-unlist(output)
    return(output)
  }
  
  ## apply select.FTD to each community in turn
  out<-apply(spmat,1,function(x) select.FTD(tdmat=tdmat,spvec=x,q=q,abund=abund))
  df.out<-data.frame(t(out))
  rownames(df.out)<-rownames(spmat)
  ## warning for zero-species communities
  if(sum(df.out$nsp==0)>0){
    warning("at least one community has no species")
  }
  
  ## calculate mean richness, dispersion, evenness, FTD
  u.M<-sum(df.out$nsp*df.out$M)/sum(df.out$nsp)
  u.nsp<-mean(df.out$nsp)
  ## to do: check if u.nsp is always calculated as arithmetic mean
  if(q==1){
    ## geometric mean -- limit of generalized mean as q->1
    u.qDT<-prod(df.out$qDT)^(1/n.comm)
  } else {
    ## generalized mean with m=1-q
    u.qDT<-(sum(df.out$qDT^(1-q))/n.comm)^1/(1-q)
  }
  u.M.prime<-u.M*u.nsp/(u.nsp-1)
  
  ## calculate mean FTD and evenness
  u.qDTM<-1+u.qDT*u.M
  u.Et<-u.qDT/u.nsp
  
  ## list more things
  list(com.FTD=df.out,u.nsp=u.nsp,u.M=u.M,u.M.prime=u.M.prime,u.Et=u.Et,u.qDT=u.qDT,u.qDTM=u.qDTM)
}

## to consider:
## change to output standardized M rather than normal
## structured vs. unstructured gamma-diversity