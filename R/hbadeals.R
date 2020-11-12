library(coda)
library(limma)
library(rstan)
library(bayestestR)

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



hbadeals.heirarchy=function(countsData,labels,n.cores=getOption("mc.cores", 2L),isoform.level=FALSE,mcmc.iter=3000,mcmc.warmup=4000,lib.size=NULL,mtc=FALSE)
{
  if (sum(labels!=labels[order(labels)])>0) {print('Error: The samples are not ordered!');return(NULL)}
  labels=factor(labels)
  design=model.matrix(~group,data.frame(group=labels))
  summed.counts=do.call(rbind,lapply(split(countsData[,3:ncol(countsData)],countsData[,1]),function(m){colSums(m)}))
  summed.counts=summed.counts[order(match(rownames(summed.counts),as.character(countsData[,1]))),]

  iso.data=getvar(countsData[3:ncol(countsData)],design,lib.size)
  gene.data=getvar(summed.counts,design,lib.size)
  
  theta.vals=theta.heirarchical(countsData,labels,100,lib.size,n.cores)
  
  theta_a=theta.vals[[1]]
  
  theta_b=theta.vals[[2]]

  modelString = "
  data {
    int<lower=0> Nisoforms;
    int<lower=0> Ncondition1;
    int<lower=0> Ncondition2;
    real mean_controls;
    vector[Ncondition1] sd_exp_controls;
    vector[Ncondition2] sd_exp_cases;
    vector[Nisoforms] sd_iso_cases[Ncondition2];
    vector[Nisoforms] sd_iso_controls[Ncondition1];
    vector[Nisoforms] counts_cases[Ncondition2];
    vector[Nisoforms] counts_controls[Ncondition1];
    }

  parameters {
    simplex[Nisoforms] frac;
    vector[Ncondition1] expression_controls;
    vector[Ncondition2] expression_cases;
    real beta;
    real intercept;
    simplex[Nisoforms] alpha;
    }

  model {"
  
  if (mtc)
  {
    
    modelString=paste0(modelString,  "target+=log_sum_exp(log(",theta_b,")+normal_lpdf(beta|0,5),log(",1-theta_b,")+normal_lpdf(beta|0,0.04));\n
                       
                                     target+=log_sum_exp(log(",1-theta_a,")+dirichlet_lpdf(alpha|rep_vector(100,Nisoforms)),
     log(",theta_a,")+dirichlet_lpdf(alpha|rep_vector(1,Nisoforms)));")
    
    
    
    
  }else{
    
    modelString=paste0(modelString, "alpha ~ dirichlet(rep_vector(1.0,Nisoforms));\nbeta ~ normal(0,5);\n")
  }
  
  modelString=paste0(modelString,"

      frac ~ dirichlet(rep_vector(1.0,Nisoforms));

      intercept ~ normal(mean_controls,5);

      expression_controls ~ normal(intercept,sd_exp_controls);

      expression_cases ~ normal(intercept+beta,sd_exp_cases);

      for ( j in 1:Nisoforms )
      {
          counts_controls[,j] ~ normal(log2(frac[j])+expression_controls,sd_iso_controls[,j]);

          counts_cases[,j] ~ normal(log2((frac[j]*alpha[j])/dot_product(frac,alpha))+expression_cases,sd_iso_cases[,j]);
      }

  }
  ")

  stan.mod = stan_model( model_code=modelString )

  tab=iso.data[[1]]

  summed.counts=gene.data[[1]]

  results.tab=do.call(rbind,mclapply(1:length(unique(countsData[,1])),function(gene.number){

    gene.rows=which(countsData[,1] == unique(countsData[,1])[gene.number])
    
    num.isoforms=length(gene.rows)

    initf = function() {
      list(frac = rowSums(2^tab[gene.rows,])/sum(rowSums(2^tab[gene.rows,])),
           expression_cases=summed.counts[gene.number,labels==2],
           expression_controls=summed.counts[gene.number,labels==1],
           alpha=array(rep(1/num.isoforms,num.isoforms),dim=num.isoforms),
           beta=0,intercept=log2(mean(2^summed.counts[gene.number,labels==1]-0.5)+0.5))
    }

    res=matrix(ncol=4,nrow=0)

    colnames(res)=c('Gene','Isoform','ExplogFC/FC','P')

    dataList = list(
      counts_cases = t(tab[gene.rows,labels==2]),
      counts_controls = t(tab[gene.rows,labels==1]),
      Nisoforms = num.isoforms,
      Ncondition1 = sum(labels==1),
      Ncondition2 = sum(labels==2),
      mean_controls=log2(mean(2^summed.counts[gene.number,labels==1]-0.5)+0.5),
      sd_iso_cases=t(sqrt(iso.data[[2]][gene.rows,labels==2])),
      sd_iso_controls=t(sqrt(iso.data[[2]][gene.rows,labels==1])),
      sd_exp_cases=sqrt(gene.data[[2]][gene.number,labels==2]),
      sd_exp_controls=sqrt(gene.data[[2]][gene.number,labels==1])

    )

    stanFit = sampling( object=stan.mod , data = dataList,cores=1 ,init=initf, chains = 1,refresh=0 ,
                        iter = mcmc.warmup+mcmc.iter,warmup=mcmc.warmup , thin = 1 )

    mcmcCoda = mcmc.list( lapply( 1:ncol(stanFit) , function(x) { mcmc(as.array(stanFit)[,x,]) } ) )

    beta.col=which(colnames(mcmcCoda[[1]])=='beta')

    mean.de=mean(mcmcCoda[[1]][,beta.col])

    expression.p=p_rope(as.numeric(mcmcCoda[[1]][,beta.col]) ,range = c(-0.1, 0.1))$p_ROPE
  
    gene.name=as.character(unique(countsData[,1])[gene.number])

    res=rbind(res,c(gene.name,'Expression',mean.de,expression.p))

    frac=rep(0,num.isoforms)

    alpha=rep(0,num.isoforms)

    alpha.p=rep(0,num.isoforms)

    for (i in (1:num.isoforms))
    {

      next.frac.col=which(colnames(mcmcCoda[[1]])==paste0('frac[',i,']'))

      frac[i]=mean(mcmcCoda[[1]][,next.frac.col])

      next.alpha.col=which(colnames(mcmcCoda[[1]])==paste0('alpha[',i,']'))

      alpha[i]=mean(mcmcCoda[[1]][,next.alpha.col])
    }

    if (isoform.level){
      uni.val=sum(frac*alpha)
    }else{
      uni.val=1/num.isoforms
    }
    for (i in (1:num.isoforms))
    {

      next.alpha.col=which(colnames(mcmcCoda[[1]])==paste0('alpha[',i,']'))

      alpha.p[i]=p_rope(as.numeric(mcmcCoda[[1]][,next.alpha.col])-uni.val,range = c(-0.1*uni.val,0.1*uni.val))$p_ROPE
    }
  
    frac.2=(frac*alpha)/sum(frac*alpha)

    fc.ds=frac.2/frac

    if (isoform.level)
    {
      for (i in (1:num.isoforms))

        res=rbind(res,c(gene.name,as.character(countsData[gene.rows[i],2]),fc.ds[i],alpha.p[i]))
    }else{
        min.p=min(alpha.p)
        res=rbind(res,c(gene.name,'Splicing',fc.ds[which(alpha.p==min(alpha.p))[1]],min.p))
    }
    return(res)
  },mc.cores = n.cores))

  return(results.tab)

}

hbadeals.flat=function(countsData,labels,n.cores=getOption("mc.cores", 2L),isoform.level=FALSE,mcmc.iter=3000,mcmc.warmup=4000,lib.size=NULL,mtc=FALSE)
{
    if (sum(labels!=labels[order(labels)])>0) {print('Error: The samples are not ordered!');return(NULL)}
    labels=factor(labels)
    design=model.matrix(~group,data.frame(group=labels))
    summed.counts=do.call(rbind,lapply(split(countsData[,3:ncol(countsData)],countsData[,1]),function(m){colSums(m)}))
    summed.counts=summed.counts[order(match(rownames(summed.counts),as.character(countsData[,1]))),]
    
    iso.data=getvar(countsData[3:ncol(countsData)],design,lib.size)
    gene.data=getvar(summed.counts,design,lib.size)
 
    theta.vals=theta.flat(countsData,labels,100,lib.size,n.cores)
    
    theta_a=theta.vals[[1]]
    
    theta_b=theta.vals[[2]]
       
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
    
    model {"

    if (mtc)
    {
  
      modelString=paste0(modelString,  "target+=log_sum_exp(log(",theta_b,")+normal_lpdf(beta|0,5),log(",1-theta_b,")+normal_lpdf(beta|0,0.04));\n
                       
                     target+=log_sum_exp(log(",1-theta_a,")+dirichlet_lpdf(alpha|rep_vector(100,Nisoforms)),
     log(",theta_a,")+dirichlet_lpdf(alpha|rep_vector(1,Nisoforms)));")
      
    }else{
      
      modelString=paste0(modelString, "alpha ~ dirichlet(rep_vector(1.0,Nisoforms));\nbeta ~ normal(0,5);\n")
    }
    
     modelString=paste0(modelString, "
        
        frac ~ dirichlet(rep_vector(1.0,Nisoforms));
      
        intercept ~ normal(mean_controls,5);
        
        for ( j in 1:Nisoforms )
        {
            counts_controls[,j] ~ normal(log2(frac[j])+intercept,sd_iso_controls[,j]);
            
            counts_cases[,j] ~ normal(log2((frac[j]*alpha[j])/dot_product(frac,alpha))+intercept+beta,sd_iso_cases[,j]);
        }
        
    }
    ")
    
    stan.mod = stan_model( model_code=modelString )
    
    tab=iso.data[[1]]
    
    summed.counts=gene.data[[1]]
    
    results.tab=do.call(rbind,mclapply(1:length(unique(countsData[,1])),function(gene.number){
        
        gene.rows=which(countsData[,1] == unique(countsData[,1])[gene.number])
        
        num.isoforms=length(gene.rows)
        
        initf = function() {
            list(frac = rowSums(2^tab[gene.rows,])/sum(rowSums(2^tab[gene.rows,])),
            alpha=array(rep(1/num.isoforms,num.isoforms),dim=num.isoforms),
            beta=0,intercept=log2(mean(2^summed.counts[gene.number,labels==1]-0.5)+0.5))
        }
        
        res=matrix(ncol=4,nrow=0)
        
        colnames(res)=c('Gene','Isoform','ExplogFC/FC','P')
        
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
        
        mcmcCoda = mcmc.list( lapply( 1:ncol(stanFit) , function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
        
        beta.col=which(colnames(mcmcCoda[[1]])=='beta')
        
        mean.de=mean(mcmcCoda[[1]][,beta.col])
        
        expression.p=p_rope(as.numeric(mcmcCoda[[1]][,beta.col]), range = c(-0.1, 0.1))$p_ROPE
        
        gene.name=as.character(unique(countsData[,1])[gene.number])
        
        res=rbind(res,c(gene.name,'Expression',mean.de,expression.p))
        
        frac=rep(0,num.isoforms)
        
        alpha=rep(0,num.isoforms)
        
        alpha.p=rep(0,num.isoforms)
        
        for (i in (1:num.isoforms))
        {
            
            next.frac.col=which(colnames(mcmcCoda[[1]])==paste0('frac[',i,']'))
            
            frac[i]=mean(mcmcCoda[[1]][,next.frac.col])
            
            next.alpha.col=which(colnames(mcmcCoda[[1]])==paste0('alpha[',i,']'))
            
            alpha[i]=mean(mcmcCoda[[1]][,next.alpha.col])
        }
        
        if (isoform.level){
            uni.val=sum(frac*alpha)
        }else{
            uni.val=1/num.isoforms
        }
        for (i in (1:num.isoforms))
        {
            
            next.alpha.col=which(colnames(mcmcCoda[[1]])==paste0('alpha[',i,']'))
            
            alpha.p[i]=p_rope(as.numeric(mcmcCoda[[1]][,next.alpha.col])-uni.val,range = c(-0.1*uni.val,0.1*uni.val))$p_ROPE
            
        }
        
        frac.2=(frac*alpha)/sum(frac*alpha)
        
        fc.ds=frac.2/frac
        
        if (isoform.level)
        {
            for (i in (1:num.isoforms))
            
            res=rbind(res,c(gene.name,as.character(countsData[gene.rows[i],2]),fc.ds[i],alpha.p[i]))
        }else{
            min.p=min(alpha.p)
            res=rbind(res,c(gene.name,'Splicing',fc.ds[which(alpha.p==min(alpha.p))[1]],min.p))
        }
        return(res)
    },mc.cores = n.cores))
    
    return(results.tab)
    
}

hbadeals=function(countsData,labels,n.cores=getOption("mc.cores", 2L),isoform.level=FALSE,mcmc.iter=3000,mcmc.warmup=4000,hierarchy='auto',lib.size=NULL,mtc=FALSE,theta_a=0.5,theta_b=0.5)
{
        use.heirarchical=TRUE

        if (hierarchy=='auto')
        {
            use.heirarchical=choose.model(countsData=countsData,labels=labels,n.cores=n.cores,num.genes=100, mcmc.iter=mcmc.iter, mcmc.warmup=mcmc.warmup,lib.size=lib.size)
            
        }else if (hierarchy=='no')
        {
           use.heirarchical=FALSE
            
        }
        
    if (use.heirarchical)
    
        return (hbadeals.heirarchy(countsData = countsData,labels = labels,n.cores=n.cores,isoform.level = isoform.level,mcmc.iter=mcmc.iter,mcmc.warmup=mcmc.warmup,lib.size=lib.size,mtc=mtc))
    
    return (hbadeals.flat(countsData = countsData,labels = labels,n.cores=n.cores,isoform.level = isoform.level,mcmc.iter=mcmc.iter,mcmc.warmup=mcmc.warmup,lib.size=lib.size,mtc=mtc))
    

}

