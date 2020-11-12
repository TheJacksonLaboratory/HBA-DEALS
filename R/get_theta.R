

theta.flat=function(countsData,labels,opt.iter=30,lib.size=NULL,n.cores=getOption("mc.cores", 2L))
{
  if (sum(labels!=labels[order(labels)])>0) {print('Error: The samples are not ordered!');return(NULL)}
  labels=factor(labels)
  design=model.matrix(~group,data.frame(group=labels))
  summed.counts=do.call(rbind,lapply(split(countsData[,3:ncol(countsData)],countsData[,1]),function(m){colSums(m)}))
  summed.counts=summed.counts[order(match(rownames(summed.counts),as.character(countsData[,1]))),]
  
  iso.data=getvar(countsData[3:ncol(countsData)],design,lib.size)
  gene.data=getvar(summed.counts,design,lib.size)
  
  modelString="data {
  int<lower=0> Ngenes;
  int<lower=0> Ntranscripts;
  int<lower=0> Ncondition1;
  int<lower=0> Ncondition2;
  int<lower=0> Nisoforms[Ngenes];
  real mean_controls[Ngenes];
  matrix[Ncondition2,Ntranscripts] sd_iso_cases;
  matrix[Ncondition1,Ntranscripts] sd_iso_controls;
  matrix[Ncondition2,Ntranscripts] counts_cases;
  matrix[Ncondition1,Ntranscripts] counts_controls;
  int<lower=0> num_opt_genes;
  int run_ids[num_opt_genes];
  
} 
  parameters { 
  vector[Ngenes] beta;
  vector[Ngenes] intercept;
  vector<lower=0,upper=1>[Ntranscripts] alpha;
  vector<lower=0,upper=1>[Ntranscripts] frac;
  real<lower=0,upper=1> theta_a;
  real<lower=0,upper=1> theta_b;
  
  }
  model {
  int first_i=1;
  int last_i=0;
  real sum_a;
  real sum_f;
  
  for (genes_itr in run_ids[1:num_opt_genes])
  {
  last_i=first_i+Nisoforms[genes_itr]-1;
  
  sum_a=sum(alpha[first_i:last_i]);
  
  sum_f=sum(frac[first_i:last_i]);
  
  target+=dirichlet_lpdf(frac[first_i:last_i]/sum_f|rep_vector(1,Nisoforms[genes_itr]));
  
  target+=log_sum_exp(log(theta_b)+normal_lpdf(beta[genes_itr]|0,5),log(1-theta_b)+normal_lpdf(beta[genes_itr]|0,0.04));
  
  
  target+=log_sum_exp(log(1-theta_a)+dirichlet_lpdf(alpha[first_i:last_i]/sum_a|rep_vector(100,Nisoforms[genes_itr])),
  log(theta_a)+dirichlet_lpdf(alpha[first_i:last_i]/sum_a|rep_vector(1,Nisoforms[genes_itr])));
  
  
  intercept[genes_itr] ~ normal(mean_controls[genes_itr],5);
  
  for (j in first_i:last_i)
  {
  target+=normal_lpdf(counts_controls[,j] | log2(frac[j]/sum_f)+intercept[genes_itr],sd_iso_controls[,j]);
  
  target+=normal_lpdf(counts_cases[,j] | log2((frac[j]/sum_f*alpha[j]/sum_a)/dot_product(frac[first_i:last_i]/sum_f,alpha[first_i:last_i]/sum_a))+intercept[genes_itr]+beta[genes_itr],sd_iso_cases[,j]);
  }
  
  first_i=last_i+1;
  }
  }"
  
  stan.mod = stan_model( model_code=modelString )
  
  tab=iso.data[[1]]
  
  summed.counts=gene.data[[1]]
  
  isoform.count=table(countsData[,1])
  
  isoform.count=isoform.count[match(unique(countsData[,1]),names(isoform.count))]
  
  isoform.count=as.numeric(isoform.count)
  
  split.factor=factor(as.character(countsData[,1]),levels=unique(as.character(countsData[,1])))
  
  init_l=list(theta_a=0.05,
              theta_b=0.05,
              frac=array(do.call(c,lapply(lapply(split(tab,split.factor),matrix,ncol=ncol(tab)),function(m)rowSums(2^m)/sum(rowSums(2^m)))),dim=nrow(countsData)),
              alpha=array(1/unlist(lapply(isoform.count,function(x)rep(x,x))),dim=nrow(countsData)),
              beta=rep(0,nrow(summed.counts)),
              intercept=log2(rowMeans(2^summed.counts[,labels==1]-0.5)+0.5)
  )
  
  gene.names=unique(countsData[,1])
  
  g.group=1:nrow(summed.counts)
  
  dataList=list(Ncondition1 = sum(labels==1),
                Ncondition2 = sum(labels==2),
                counts_cases = t(tab[,labels==2]),
                counts_controls = t(tab[,labels==1]),
                mean_controls=log2(rowMeans(2^summed.counts[,labels==1]-0.5)+0.5),
                sd_iso_cases=t(sqrt(iso.data[[2]][,labels==2])),
                sd_iso_controls=t(sqrt(iso.data[[2]][,labels==1])),
                Ngenes=nrow(summed.counts),
                Ntranscripts=nrow(countsData),
                Nisoforms=isoform.count,
                run_ids=g.group,
                num_opt_genes=length(g.group)
  )
  
  stanFit = optimizing( object=stan.mod , data = dataList,init=init_l, iter = opt.iter)
  
  theta_a=stanFit$par[grepl('theta_a',names(stanFit$par))]
  
  theta_b=stanFit$par[grepl('theta_b',names(stanFit$par))]
  
  return(list(theta_a,theta_b))
  
  }






theta.heirarchical=function(countsData,labels,opt.iter=30,lib.size=NULL,n.cores=getOption("mc.cores", 2L))
{
  if (sum(labels!=labels[order(labels)])>0) {print('Error: The samples are not ordered!');return(NULL)}
  labels=factor(labels)
  design=model.matrix(~group,data.frame(group=labels))
  summed.counts=do.call(rbind,lapply(split(countsData[,3:ncol(countsData)],countsData[,1]),function(m){colSums(m)}))
  summed.counts=summed.counts[order(match(rownames(summed.counts),as.character(countsData[,1]))),]
  
  iso.data=getvar(countsData[3:ncol(countsData)],design,lib.size)
  gene.data=getvar(summed.counts,design,lib.size)
  
  modelString="data {
  int<lower=0> Ngenes;
  int<lower=0> Ntranscripts;
  int<lower=0> Ncondition1;
  int<lower=0> Ncondition2;
  int<lower=0> Nisoforms[Ngenes];
  real mean_controls[Ngenes];
  matrix[Ngenes,Ncondition2] sd_exp_cases;
  matrix[Ngenes,Ncondition1] sd_exp_controls;
  matrix[Ncondition2,Ntranscripts] sd_iso_cases;
  matrix[Ncondition1,Ntranscripts] sd_iso_controls;
  matrix[Ncondition2,Ntranscripts] counts_cases;
  matrix[Ncondition1,Ntranscripts] counts_controls;
  int<lower=0> num_opt_genes;
  int run_ids[num_opt_genes];
  
} 
  parameters { 
  vector[Ngenes] beta;
  vector[Ngenes] intercept;
  vector<lower=0,upper=1>[Ntranscripts] alpha;
  vector<lower=0,upper=1>[Ntranscripts] frac;
  matrix[Ngenes,Ncondition1] expression_controls;
  matrix[Ngenes,Ncondition2] expression_cases;
  real<lower=0,upper=1> theta_a;
  real<lower=0,upper=1> theta_b;
  
  }
  model {
  int first_i=1;
  int last_i=0;
  real sum_a;
  real sum_f;
  
  
  for (genes_itr in run_ids[1:num_opt_genes])
  {
  last_i=first_i+Nisoforms[genes_itr]-1;
  
  sum_a=sum(alpha[first_i:last_i]);
  
  sum_f=sum(frac[first_i:last_i]);
  
  target+=dirichlet_lpdf(frac[first_i:last_i]/sum_f|rep_vector(1,Nisoforms[genes_itr]));
  
  target+=log_sum_exp(log(theta_b)+normal_lpdf(beta[genes_itr]|0,5),log(1-theta_b)+normal_lpdf(beta[genes_itr]|0,0.04));
  
  target+=log_sum_exp(log(1-theta_a)+dirichlet_lpdf(alpha[first_i:last_i]/sum_a|rep_vector(100,Nisoforms[genes_itr])),
     log(theta_a)+dirichlet_lpdf(alpha[first_i:last_i]/sum_a|rep_vector(1,Nisoforms[genes_itr])));
  
  
  intercept[genes_itr] ~ normal(mean_controls[genes_itr],5);
  
  expression_controls[genes_itr,] ~ normal(intercept[genes_itr],sd_exp_controls[genes_itr,]);
  
  expression_cases[genes_itr,] ~ normal(intercept[genes_itr]+beta[genes_itr],sd_exp_cases[genes_itr,]);
  
  for (j in first_i:last_i)
  {
  target+=normal_lpdf(counts_controls[,j] | log2(frac[j]/sum_f)+expression_controls[genes_itr,],sd_iso_controls[,j]);
  
  target+=normal_lpdf(counts_cases[,j] | log2((frac[j]/sum_f*alpha[j]/sum_a)/dot_product(frac[first_i:last_i]/sum_f,alpha[first_i:last_i]/sum_a))+expression_cases[genes_itr,],sd_iso_cases[,j]);
  }
  
  first_i=last_i+1;
  }
  }"
  
  stan.mod = stan_model( model_code=modelString )
  
  tab=iso.data[[1]]
  
  summed.counts=gene.data[[1]]
  
  isoform.count=table(countsData[,1])
  
  isoform.count=isoform.count[match(unique(countsData[,1]),names(isoform.count))]
  
  isoform.count=as.numeric(isoform.count)
  
  split.factor=factor(as.character(countsData[,1]),levels=unique(as.character(countsData[,1])))
  
  init_l=list(theta_a=0.05,
              theta_b=0.05,
              frac=array(do.call(c,lapply(lapply(split(tab,split.factor),matrix,ncol=ncol(tab)),function(m)rowSums(2^m)/sum(rowSums(2^m)))),dim=nrow(countsData)),
              alpha=array(1/unlist(lapply(isoform.count,function(x)rep(x,x))),dim=nrow(countsData)),
              beta=rep(0,nrow(summed.counts)),
              intercept=log2(rowMeans(2^summed.counts[,labels==1]-0.5)+0.5),
              expression_cases=summed.counts[,labels==2],
              expression_controls=summed.counts[,labels==1]
  )
  
  gene.names=unique(countsData[,1])
  
  g.group=1:nrow(summed.counts)
  
  dataList=list(Ncondition1 = sum(labels==1),
                Ncondition2 = sum(labels==2),
                counts_cases = t(tab[,labels==2]),
                counts_controls = t(tab[,labels==1]),
                mean_controls=log2(rowMeans(2^summed.counts[,labels==1]-0.5)+0.5),
                sd_iso_cases=t(sqrt(iso.data[[2]][,labels==2])),
                sd_iso_controls=t(sqrt(iso.data[[2]][,labels==1])),
                sd_exp_cases=sqrt(gene.data[[2]][,labels==2]),
                sd_exp_controls=sqrt(gene.data[[2]][,labels==1]),
                Ngenes=nrow(summed.counts),
                Ntranscripts=nrow(countsData),
                Nisoforms=isoform.count,
                run_ids=g.group,
                num_opt_genes=length(g.group)
                
  )
  
  stanFit = optimizing( object=stan.mod , data = dataList,init=init_l, iter = opt.iter)
  
  theta_a=stanFit$par[grepl('theta_a',names(stanFit$par))]
  
  theta_b=stanFit$par[grepl('theta_b',names(stanFit$par))]
  
  return(list(theta_a,theta_b))
  
  }


