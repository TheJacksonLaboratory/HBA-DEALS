library(bridgesampling)

choose.model=function(countsData,labels,n.cores=getOption("mc.cores", 2L),num.genes=10, mcmc.iter = 3000, mcmc.warmup = 4000,lib.size=NULL)
{
  
  sample.genes=sample(length(unique(countsData[,1])),num.genes)

  flat.mls=ml.flat(countsData,labels,n.cores,sample.genes,mcmc.iter,mcmc.warmup,lib.size)
  
  heirarchical.mls=ml.heirarchical(countsData,labels,n.cores,sample.genes,mcmc.iter,mcmc.warmup,lib.size)
  
  bf.vals=mclapply(1:num.genes,function(i){bf(flat.mls[[i]],heirarchical.mls[[i]])$bf},mc.cores = n.cores)
  
  if (median(unlist(bf.vals))>=1)
    
    return(FALSE)
  
  return(TRUE)
  
}
