

ml.heirarchical=function (countsData, labels, n.cores = getOption("mc.cores", 
                                                                  2L), gene.numbers, mcmc.iter = 3000, mcmc.warmup = 4000,lib.size=NULL) {
  if (sum(labels != labels[order(labels)]) > 0) {
    print("Error: The samples are not ordered!")
    return(NULL)
  }
  labels = factor(labels)
  design = model.matrix(~group, data.frame(group = labels))
  summed.counts = do.call(rbind, lapply(split(countsData[, 
                                                         3:ncol(countsData)], countsData[, 1]), function(m) {
                                                           colSums(m)
                                                         }))
  summed.counts = summed.counts[order(match(rownames(summed.counts), 
                                            as.character(countsData[, 1]))), ]
  iso.data = getvar(countsData[3:ncol(countsData)], design,lib.size)
  gene.data = getvar(summed.counts,design,lib.size)
  modelString = "\n  data {\n    int<lower=0> Nisoforms;\n    int<lower=0> Ncondition1;\n    int<lower=0> Ncondition2;\n    real mean_controls;\n    vector[Ncondition1] sd_exp_controls;\n    vector[Ncondition2] sd_exp_cases;\n    vector[Nisoforms] sd_iso_cases[Ncondition2];\n    vector[Nisoforms] sd_iso_controls[Ncondition1];\n    vector[Nisoforms] counts_cases[Ncondition2];\n    vector[Nisoforms] counts_controls[Ncondition1];\n    }\n\n  parameters {\n    simplex[Nisoforms] frac;\n    vector[Ncondition1] expression_controls;\n    vector[Ncondition2] expression_cases;\n    real beta;\n    real intercept;\n    simplex[Nisoforms] alpha;\n    }\n\n  model {\n\n      target+=  dirichlet_lpdf(frac|rep_vector(1.0,Nisoforms));\n\n      target+= dirichlet_lpdf(alpha|rep_vector(1.0,Nisoforms));\n\n      target+=normal_lpdf(beta |0,5);\n\n      target+= normal_lpdf(intercept|mean_controls,5);\n\n      target+=normal_lpdf(expression_controls | intercept,sd_exp_controls);\n\n      target+=normal_lpdf(expression_cases | intercept+beta,sd_exp_cases);\n\n      for ( j in 1:Nisoforms )\n      {\n          target+= normal_lpdf(counts_controls[,j]|log2(frac[j])+expression_controls,sd_iso_controls[,j]);\n\n          target+=normal_lpdf(counts_cases[,j] |log2((frac[j]*alpha[j])/dot_product(frac,alpha))+expression_cases,sd_iso_cases[,j]);\n      }\n\n  }\n  "
  stan.mod = stan_model(model_code = modelString)
  tab = iso.data[[1]]
  summed.counts = gene.data[[1]]
  
  
  res=mclapply(gene.numbers,function(gene.number)
  {
  gene.rows = which(countsData[, 1] == unique(countsData[,1])[gene.number])
  num.isoforms = length(gene.rows)
  initf = function() {
   list(frac = rowSums(2^tab[gene.rows, ])/sum(rowSums(2^tab[gene.rows,])), expression_cases = summed.counts[gene.number, labels == 2], 
        expression_controls = summed.counts[gene.number, labels == 1], alpha = array(rep(1/num.isoforms, num.isoforms), dim = num.isoforms), beta = 0, 
         intercept = log2(mean(2^summed.counts[gene.number, labels == 1] - 0.5) + 0.5))
    }
     dataList = list(counts_cases = t(tab[gene.rows, labels ==  2]), 
                     counts_controls = t(tab[gene.rows, labels == 1]), 
                     Nisoforms = num.isoforms, Ncondition1 = sum(labels == 1), 
                     Ncondition2 = sum(labels == 2), 
                     mean_controls = log2(mean(2^summed.counts[gene.number,labels == 1] - 0.5) + 0.5), 
                     sd_iso_cases = t(sqrt(iso.data[[2]][gene.rows,labels == 2])), 
                     sd_iso_controls = t(sqrt(iso.data[[2]][gene.rows,  labels == 1])), 
                     sd_exp_cases = sqrt(gene.data[[2]][gene.number, labels == 2]),
                     sd_exp_controls = sqrt(gene.data[[2]][gene.number,labels == 1]))
      stanFit = sampling(object = stan.mod, data = dataList,  cores = 1, init = initf, chains = 1, refresh = 0, 
                          iter = mcmc.warmup + mcmc.iter, warmup = mcmc.warmup, thin = 1)

      return(bridge_sampler(stanFit,silent=TRUE))
    },mc.cores = n.cores)
  return(res)                                                                   
}