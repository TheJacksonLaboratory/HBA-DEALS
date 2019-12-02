library(coda)
library(limma)
library(rstan)

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



hbadeals=function(countsData,labels,n.cores=getOption("mc.cores", 2L),isoform.level=FALSE,mcmc.iter=3000,mcmc.warmup=4000)
{
  if (sum(labels!=labels[order(labels)])>0) {print('Error: The samples are not ordered!');return(NULL)}
  labels=factor(labels)
  design=model.matrix(~group,data.frame(group=labels))
  summed.counts=do.call(rbind,lapply(split(countsData[,3:ncol(countsData)],countsData[,1]),function(m){colSums(m)}))
  summed.counts=summed.counts[order(match(rownames(summed.counts),as.character(countsData[,1]))),]

  iso.data=getvar(countsData[3:ncol(countsData)],design)
  gene.data=getvar(summed.counts,lib.size=colSums(summed.counts),design)

  modelString = "
  data {
    int<lower=0> Nisoforms;
    int<lower=0> Ncondition1;
    int<lower=0> Ncondition2;
    real mean_controls;
    real<lower=0> sd_exp[Ncondition1+Ncondition2];
    real<lower=0> sd_iso[Nisoforms,Ncondition1+Ncondition2];
    real counts[Nisoforms,Ncondition1+Ncondition2];
  }

  parameters {
    simplex[Nisoforms] frac;
    real expression[Ncondition1+Ncondition2];
    real beta;
    real intercept;
    simplex[Nisoforms] alpha;
  }

  model {

    vector[Nisoforms] dweights;

    real Z=0;

    for ( i in 1:Nisoforms )
    {
       dweights[i]=1;

       Z=Z+frac[i]*alpha[i];
    }

    frac ~ dirichlet(dweights);

    alpha ~ dirichlet(dweights);

    beta ~ normal(0,5);

    intercept ~ normal(mean_controls,5);

    for ( i in 1:(Ncondition1+Ncondition2) )
      if (i<=Ncondition1)
          expression[i] ~ normal(intercept,sd_exp[i]);
      else
          expression[i] ~ normal(intercept+beta,sd_exp[i]);

    for ( j in 1:Nisoforms )
    {
        for ( i in 1:Ncondition1 )

           counts[j,i] ~ normal(log2(frac[j])+expression[i],sd_iso[j,i]);

        for ( i in (Ncondition1+1):(Ncondition1+Ncondition2) )

           counts[j,i] ~ normal(log2((frac[j]*alpha[j])/Z)+expression[i],sd_iso[j,i]);
    }
  }
  "

  stan.mod = stan_model( model_code=modelString )

  tab=iso.data[[1]]

  summed.counts=gene.data[[1]]

  results.tab=do.call(rbind,mclapply(1:length(unique(countsData[,1])),function(gene.number){

    gene.rows=which(countsData[,1] == unique(countsData[,1])[gene.number])

    initf = function() {
      list(frac = rowSums(2^tab[gene.rows,])/sum(rowSums(2^tab[gene.rows,])),
           expression=summed.counts[gene.number,],
           alpha=array(rep(1/length(gene.rows),length(gene.rows)),dim=length(gene.rows)),
           beta=0,intercept=log2(mean(2^summed.counts[gene.number,labels==1]-0.5)+0.5))
    }

    res=matrix(ncol=4,nrow=0)

    colnames(res)=c('Gene','Isoform','ExplogFC/FC','P')

    dataList = list(
      counts = tab[gene.rows,] ,
      Nisoforms = length(gene.rows),
      Ncondition1 = sum(labels==1),
      Ncondition2 = sum(labels==2),
      mean_controls=log2(mean(2^summed.counts[gene.number,labels==1]-0.5)+0.5),
      sd_exp=sqrt(gene.data[[2]][gene.number,]),
      sd_iso=sqrt(iso.data[[2]][gene.rows,])
    )

    stanFit = sampling( object=stan.mod , data = dataList,cores=1 ,init=initf, chains = 1,refresh=0 ,iter = mcmc.warmup+mcmc.iter,warmup=mcmc.warmup , thin = 1 )

    mcmcCoda = mcmc.list( lapply( 1:ncol(stanFit) , function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

    beta.col=which(colnames(mcmcCoda[[1]])=='beta')

    mean.de=mean(mcmcCoda[[1]][,beta.col])

    if (mean.de>0)
    {
      expression.p=sum(mcmcCoda[[1]][,beta.col]<0)/mcmc.iter
    }else{

      expression.p=sum(mcmcCoda[[1]][,beta.col]>0)/mcmc.iter
    }

    fc.de=mean.de

    res=rbind(res,c(as.character(unique(countsData[,1])[gene.number]),'Expression',fc.de,expression.p))

    frac=rep(0,length(gene.rows))

    alpha=rep(0,length(gene.rows))

    alpha.p=rep(0,length(gene.rows))

    for (i in (1:length(gene.rows)))
    {

      next.frac.col=which(colnames(mcmcCoda[[1]])==paste0('frac[',i,']'))

      frac[i]=mean(mcmcCoda[[1]][,next.frac.col])

      next.alpha.col=which(colnames(mcmcCoda[[1]])==paste0('alpha[',i,']'))

      alpha[i]=mean(mcmcCoda[[1]][,next.alpha.col])
    }

    if (isoform.level){
      uni.val=sum(frac*alpha)
    }else{
      uni.val=1/length(gene.rows)
    }
    for (i in (1:length(gene.rows)))
    {

      next.alpha.col=which(colnames(mcmcCoda[[1]])==paste0('alpha[',i,']'))

      if (alpha[i]>uni.val)
      {
        alpha.p[i]=sum(mcmcCoda[[1]][,next.alpha.col]<uni.val)/mcmc.iter
      }else{
        alpha.p[i]=sum(mcmcCoda[[1]][,next.alpha.col]>uni.val)/mcmc.iter
      }

    }

    frac.2=(frac*alpha)/sum(frac*alpha)

    fc.ds=frac.2/frac

    if (isoform.level)
    {
      for (i in (1:length(gene.rows)))

        res=rbind(res,c(as.character(unique(countsData[,1])[gene.number]),as.character(countsData[gene.rows[i],2]),fc.ds[i],alpha.p[i]))
    }else{
        res=rbind(res,c(as.character(unique(countsData[,1])[gene.number]),'Splicing',fc.ds[which(alpha.p==min(alpha.p))[1]],min(alpha.p)))
    }
    return(res)
  },mc.cores = n.cores))

  return(results.tab)

}

