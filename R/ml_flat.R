
library(coda)
library(limma)
library(rstan)
library(bridgesampling)

getvar=function (counts, design, lib.size = NULL, normalize.method = "none")
{
  
  counts = as.matrix(counts)
  n = nrow(counts)
  if (n < 2L)
    stop("Need at least two genes")
  if (is.null(design))
    stop("Design matrix was not provided")
  if (is.null(lib.size))
    lib.size = colSums(counts)
  y = t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y = normalizeBetweenArrays(y, method = normalize.method)
  fit = lmFit(y, design)
  if (is.null(fit$Amean))
    fit$Amean = rowMeans(y, na.rm = TRUE)
  sx = fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy = sqrt(fit$sigma)
  allzero = rowSums(counts) == 0
  if (any(allzero)) {
    sx = sx[!allzero]
    sy = sy[!allzero]
  }
  l = lowess(sx, sy, f = 0.5)
  f = approxfun(l, rule = 2)
  if (fit$rank < ncol(design)) {
    j = fit$pivot[1:fit$rank]
    fitted.values = fit$coef[, j, drop = FALSE] %*% t(fit$design[, j, drop = FALSE])
  }else {
    fitted.values = fit$coef %*% t(fit$design)
  }
  fitted.cpm = 2^fitted.values
  fitted.count = 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount = log2(fitted.count)
  w = f(fitted.logcount)^4
  dim(w) = dim(fitted.logcount)
  return(list(y,w))
}



ml.flat=function(countsData,labels,n.cores=getOption("mc.cores", 2L),gene.numbers,mcmc.iter=3000,mcmc.warmup=4000,lib.size=NULL)
{
  if (sum(labels!=labels[order(labels)])>0) {print('Error: The samples are not ordered!');return(NULL)}
  labels=factor(labels)
  design=model.matrix(~group,data.frame(group=labels))
  summed.counts=do.call(rbind,lapply(split(countsData[,3:ncol(countsData)],countsData[,1]),function(m){colSums(m)}))
  summed.counts=summed.counts[order(match(rownames(summed.counts),as.character(countsData[,1]))),]
  
  iso.data=getvar(countsData[3:ncol(countsData)],design,lib.size)
  gene.data=getvar(summed.counts,design,lib.size)
  
  modelString = "
  data {
  int<lower=0> Nisoforms;
  int<lower=0> Ncondition1;
  int<lower=0> Ncondition2;
  real mean_controls;
  vector[Nisoforms] sd_iso_cases[Ncondition2];
  vector[Nisoforms] sd_iso_controls[Ncondition1];
  vector[Nisoforms] counts_cases[Ncondition2];
  vector[Nisoforms] counts_controls[Ncondition1];
  }
  
  parameters {
  simplex[Nisoforms] frac;
  real beta;
  real intercept;
  simplex[Nisoforms] alpha;
  }
  
  model {
  
  target += dirichlet_lpdf(frac|rep_vector(1.0,Nisoforms));
  
  target +=dirichlet_lpdf(alpha|rep_vector(1.0,Nisoforms));
  
  target += normal_lpdf(beta|0,5);
  
  target+= normal_lpdf(intercept|mean_controls,5);
  
  for ( j in 1:Nisoforms )
  {
    target +=  normal_lpdf(counts_controls[,j] | log2(frac[j])+intercept,sd_iso_controls[,j]);
  
    target += normal_lpdf(counts_cases[,j] | log2((frac[j]*alpha[j])/dot_product(frac,alpha))+intercept+beta,sd_iso_cases[,j]);
  }
  
  }
  "
  
  stan.mod = stan_model( model_code=modelString )
  
  tab=iso.data[[1]]
  
  summed.counts=gene.data[[1]]
  
  res=mclapply(gene.numbers,function(gene.number)
  {
    gene.rows=which(countsData[,1] == unique(countsData[,1])[gene.number])
    
    num.isoforms=length(gene.rows)
    
    initf = function() {
      list(frac = rowSums(2^tab[gene.rows,])/sum(rowSums(2^tab[gene.rows,])),
           alpha=array(rep(1/num.isoforms,num.isoforms),dim=num.isoforms),
           beta=0,intercept=log2(mean(2^summed.counts[gene.number,labels==1]-0.5)+0.5))
    }
     
    dataList = list(
      counts_cases = t(tab[gene.rows,labels==2]),
      counts_controls = t(tab[gene.rows,labels==1]),
      Nisoforms = num.isoforms,
      Ncondition1 = sum(labels==1),
      Ncondition2 = sum(labels==2),
      mean_controls=log2(mean(2^summed.counts[gene.number,labels==1]-0.5)+0.5),
      sd_iso_cases=t(sqrt(iso.data[[2]][gene.rows,labels==2])),
      sd_iso_controls=t(sqrt(iso.data[[2]][gene.rows,labels==1]))
      
    )
    
    stanFit = sampling( object=stan.mod , data = dataList,cores=1 ,init=initf, chains = 1,refresh=0 ,
                        iter = mcmc.warmup+mcmc.iter,warmup=mcmc.warmup , thin = 1 )
    
    return(bridge_sampler(stanFit,silent=TRUE))
  },mc.cores = n.cores)
  
  return(res)
}

