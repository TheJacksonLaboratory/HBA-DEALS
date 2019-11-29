library(coda)
library(limma)
library(rstan)
library(matrixStats)
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


norm.to.frac=function(v){c((exp(v))/(sum(exp(v))+1),1/(sum(exp(v))+1))}

frac.to.norm=function(v){log(v[1:(length(v)-1)]/v[length(v)])}


hbadeals=function(countsData,labels,n.cores=getOption("mc.cores", 2L),isoform.level=FALSE,mcmc.iter=1000,mcmc.warmup=500)
{
  if (sum(labels!=labels[order(labels)])>0) {print('Error: The samples are not ordered!');return(NULL)}
  labels=factor(labels)
  design=model.matrix(~group,data.frame(group=labels))
  summed.counts=do.call(rbind,lapply(split(countsData[,3:ncol(countsData)],countsData[,1]),function(m){colSums(m)}))
  summed.counts=summed.counts[order(match(rownames(summed.counts),as.character(countsData[,1]))),]

  iso.data=get.var(countsData[3:ncol(countsData)],design)
  gene.data=get.var(summed.counts,lib.size=colSums(summed.counts),design)

  modelString = "
  data {
    int<lower=0> Nisoforms;
    int<lower=0> Ncondition1;
    int<lower=0> Ncondition2;
    real exp_mean_controls;
    real<lower=0> sd_exp[Ncondition1+Ncondition2];
    real<lower=0> sd_iso[Nisoforms,Ncondition1+Ncondition2];
    real counts[Nisoforms,Ncondition1+Ncondition2];
    real normed_frac_means[Nisoforms-1];
    real normed_frac_sds1[Nisoforms-1];
    real normed_frac_sds2[Nisoforms-1];
  }

  parameters {

    matrix [Nisoforms-1,Ncondition1+Ncondition2] normed_fracs;
    real expression[Ncondition1+Ncondition2];
    real beta;
    real intercept;
    real alpha[Nisoforms-1];
  }

  transformed parameters{

    simplex[Nisoforms] fracs[Ncondition1+Ncondition2];

    for (j in 1:(Ncondition1+Ncondition2))
    {
      real summed_expo=0;

      for (i in 1:(Nisoforms-1))

         summed_expo+=exp(normed_fracs[i,j]);

      for (i in 1:(Nisoforms-1))

         fracs[j,i]=(exp(normed_fracs[i,j]))/(summed_expo+1);

      fracs[j,Nisoforms]=1/(summed_expo+1);
    }

  }

  model {

    for (j in 1:(Nisoforms-1))

       alpha[j] ~ normal(0,5);

    beta ~ normal(0,5);

    for (j in 1:(Nisoforms-1))

      for ( i in 1:(Ncondition1+Ncondition2) )
      {
        if (i<=Ncondition1)

            normed_fracs[j,i] ~ normal(normed_frac_means[j],normed_frac_sds1[j]);

        else

            normed_fracs[j,i] ~ normal(normed_frac_means[j]+alpha[j],normed_frac_sds2[j]);

      }

    intercept ~ normal(exp_mean_controls,5);

    for ( i in 1:(Ncondition1+Ncondition2) )

      if (i<=Ncondition1)

         expression[i] ~ normal(intercept,sd_exp[i]);

       else

         expression[i] ~ normal(intercept+beta,sd_exp[i]);

    for ( j in 1:Nisoforms )

      for ( i in 1:(Ncondition1+Ncondition2) )

        counts[j,i] ~ normal(log2(fracs[i,j])+expression[i],sd_iso[j,i]);

  }
  "

  stan.mod = stan_model( model_code=modelString )

  tab=iso.data[[1]]

  summed.counts=gene.data[[1]]

  results.tab=do.call(rbind,mclapply(1:length(unique(countsData[,1])),function(gene.number){

    gene.rows=which(countsData[,1] == unique(countsData[,1])[gene.number])

    m.frac=t(t(2^tab[gene.rows,])/colSums(2^tab[gene.rows,]))

    m.norm=apply(m.frac,2,frac.to.norm)

    if (length(gene.rows)==2)

      m.norm=matrix(m.norm,nrow=1)

    if (length(gene.rows)>2)
    {
      normed.frac.means=array(rowMeans(m.norm[,labels==1]),dim=length(gene.rows)-1)

      normed.frac.sds1=array(rowSds(m.norm[,labels==1]),dim=length(gene.rows)-1)

      normed.frac.sds2=array(rowSds(m.norm[,labels==2]),dim=length(gene.rows)-1)
    }else{

      normed.frac.means=array(mean(m.norm[1,labels==1]),dim=1)

      normed.frac.sds1=array(sd(m.norm[1,labels==1]),dim=1)

      normed.frac.sds2=array(sd(m.norm[1,labels==2]),dim=1)

    }

    normed.frac.sds1[normed.frac.sds1==0]=10^-3

    normed.frac.sds2[normed.frac.sds2==0]=10^-3

    initf = function() {
      list(
           expression=summed.counts[gene.number,],
           intercept=log2(mean(2^summed.counts[gene.number,labels==1]-0.5)+0.5),
           normed_fracs=m.norm,
           alpha=array(rep(0,length(gene.rows)-1),dim=length(gene.rows)-1),
           beta=0
           )
    }

    res=matrix(ncol=4,nrow=0)

    colnames(res)=c('Gene','Isoform','ExplogFC/FC','P')

    dataList = list(
      counts = tab[gene.rows,] ,
      Nisoforms = length(gene.rows),
      Ncondition1 = sum(labels==1),
      Ncondition2 = sum(labels==2),
      exp_mean_controls=log2(mean(2^summed.counts[gene.number,labels==1]-0.5)+0.5),
      sd_exp=sqrt(gene.data[[2]][gene.number,]),
      sd_iso=sqrt(iso.data[[2]][gene.rows,]),
      normed_frac_means=normed.frac.means,
      normed_frac_sds1=normed.frac.sds1,
      normed_frac_sds2=normed.frac.sds2
    )

    stanFit = sampling( object=stan.mod , data = dataList,cores=1 ,init=initf, chains = 1,refresh=0 ,iter = mcmc.warmup+mcmc.iter,warmup=mcmc.warmup,
                        thin = 1 ,pars=c('fracs','beta','alpha'))

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

    mean.fracs.control=matrix(ncol=length(gene.rows),nrow=0)

    frac.cols.control=c()

    for (i in (1:sum(labels==1)))

       frac.cols.control=c(frac.cols.control,which(grepl(paste0('fracs\\[',i,',[0-9]*\\]'),colnames(mcmcCoda[[1]]))))

    alpha=matrix(ncol=length(gene.rows)-1,nrow=0)

    alpha.cols=which(grepl('alpha\\[',colnames(mcmcCoda[[1]])))

    for (i in (1:nrow(mcmcCoda[[1]])))
    {
      mean.on=unlist(rep(1:length(gene.rows),sum(labels==1)))

      mean.fracs.control=rbind(mean.fracs.control,unlist(lapply(1:(length(gene.rows)),function(i,v,l){mean(v[l==i])},mcmcCoda[[1]][i,frac.cols.control],mean.on)))

      alpha=rbind(alpha,mcmcCoda[[1]][i,alpha.cols])
    }
    #compare to alpha[i]:
    if (length(gene.rows)>2)
    {
      uni.val=log(sum(colMeans(mean.fracs.control[,1:(length(gene.rows)-1)])*exp(colMeans(alpha)))+mean(mean.fracs.control[,length(gene.rows)]))
    }else{
      uni.val=log(sum(mean(mean.fracs.control[,1:(length(gene.rows)-1)])*exp(mean(alpha)))+mean(mean.fracs.control[,length(gene.rows)]))
    }

    mean.for.last=mean(rowSums(mean.fracs.control[,1:(length(gene.rows)-1)]*exp(alpha)))  #compare to 1-mean.fracs.control[,length(gene.rows)]

    last.uni.val=mean(1-mean.fracs.control[,length(gene.rows)])

    alpha.p=rep(0,length(gene.rows))

    for (i in (1:length(gene.rows)))
    {
      if (i<length(gene.rows))
      {
        if (mean(alpha[,i])>uni.val)
        {
          alpha.p[i]=sum(alpha[,i]<uni.val)/mcmc.iter
        }else{
          alpha.p[i]=sum(alpha[,i]>uni.val)/mcmc.iter
        }
      }else{
        if (mean.for.last>last.uni.val)
        {
          alpha.p[i]=sum(rowSums(mean.fracs.control[,1:(length(gene.rows)-1)]*exp(alpha))<last.uni.val)/mcmc.iter
        }else{
          alpha.p[i]=sum(rowSums(mean.fracs.control[,1:(length(gene.rows)-1)]*exp(alpha))>last.uni.val)/mcmc.iter
        }
      }
    }

    fc.ds=norm.to.frac(frac.to.norm(colMeans(mean.fracs.control))+colMeans(alpha))/colMeans(mean.fracs.control)

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

