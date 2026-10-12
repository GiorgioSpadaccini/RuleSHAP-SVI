RubinGelm=function(models){
m=length(models)
n=length(models[[1]]$Sigma2Samples)
P=nrow(models[[1]]$BetaSamples)+1 #Also include sigma2
Between=matrix(0,ncol=P,nrow=P)
Within=matrix(0,ncol=P,nrow=P)
grand_mu=rep(0,P)

#Compute within-chain variance mat
for(j in 1:m){
   chain_mu=c(models[[j]]$BetaHat,models[[j]]$Sigma2Hat)
   grand_mu=grand_mu+chain_mu
   Within=Within+crossprod(t(rbind(models[[j]]$BetaSamples,models[[j]]$Sigma2Samples)-chain_mu))
}
grand_mu=grand_mu/m
Within=Within/m/(n-1)

#Compute between-chain variance mat
for(j in 1:m){
   chain_mu=c(models[[j]]$BetaHat,models[[j]]$Sigma2Hat)
   Between=Between+crossprod(t(chain_mu-grand_mu))
}
Between=Between*n/(m-1)

output=c((n-1+(1+1/m)*RSpectra::eigs(solve(Within,Between), k = 1, which = "LM")$values)/n,RSpectra::eigs(Within, k = 1, which = "LM")$values)
return(output)
}



RG_conv=function(models,min_length=NULL,ping=1e3){
  n=length(models[[1]]$Sigma2Samples)
  n_iter=ceiling(n/ping)
  if(is.null(min_length)){
    m=length(models)
    P=nrow(models[[1]]$BetaSamples)+1 #Also include sigma2
    min_length=ceiling(P/m)
  }
  out=matrix(ncol=2,nrow=n_iter)
  
  for(i in n_iter:1){
    chain_length=min(i*ping,n)
    
    if(chain_length < min_length){next}
    
    models_trunc=lapply(models,function(mod){list(BetaSamples=mod$BetaSamples[,1:chain_length],
                                                  BetaHat=rowMeans(mod$BetaSamples[,1:chain_length]),
                                                  Sigma2Samples=mod$Sigma2Samples[1:chain_length],
                                                  Sigma2Hat=mean(mod$Sigma2Samples[1:chain_length]))})
    out[i,]=RubinGelm(models_trunc)
    print(paste0('MPSRF after ',chain_length,' samples: ',out[i,1],'; Max. variance within chains: ',out[i,2]))
  }
  
  return(out)
}
