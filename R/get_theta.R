write.block=function (n.iso) 
{
  paste0("\n  idx= gene_to_ds[genes_itr]; \n\n  target+=dirichlet_lpdf(frac_", 
         n.iso, "[idx]|rep_vector(1,", n.iso, "));\n  \n  target+=dirichlet_lpdf_m(alpha_", 
         n.iso, "[idx],theta_a);\n\n  l=1;\n\n  for (j in first_i:last_i)\n  {\n\n  target+=normal_lpdf(counts_controls[,j] | log2(frac_", 
         n.iso, "[idx][l])+intercept[genes_itr],sd_iso_controls[,j]);\n  \n  target+=normal_lpdf(counts_cases[,j] | log2((frac_", 
         n.iso, "[idx][l]*alpha_", n.iso, "[idx][l])/dot_product(frac_", 
         n.iso, "[idx],alpha_", n.iso, "[idx]))+intercept[genes_itr]+beta[genes_itr],sd_iso_cases[,j]);\n\n  l=l+1;\n  }\n  \n\n")
}


theta.flat=function (countsData, labels, opt.iter = 30, lib.size = NULL, 
                     n.cores = getOption("mc.cores", 2L)) 
{
  if (sum(labels != labels[order(labels)]) > 0) {
    print("Error: The samples are not ordered!")
    return(NULL)
  }
  labels = factor(labels)
  
  
  functions.string <- "
  
  functions{
  
  
  real normal_lpdf_m(real y,real theta){
  
  real s=normal_lpdf(y|0,0.04);
  
  real l=normal_lpdf(y|0,5);
  
  if (l>s)
  
  return log(theta)+l;
  
  return log(1-theta)+s;
  
  }
  
  
  real dirichlet_lpdf_m(vector y,real theta){
  
  real s=dirichlet_lpdf(y|rep_vector(50.0,num_elements(y)));
  
  real l=dirichlet_lpdf(y|rep_vector(1.0,num_elements(y)));
  
  if (l>s)
  
  return log(theta)+l;
  
  return log(1-theta)+s;
  
  }
  
  
  
  }\n"
  
  
  design = model.matrix(~group, data.frame(group = labels))
  summed.counts = do.call(rbind, lapply(split(countsData[, 
                                                         3:ncol(countsData)], countsData[, 1]), function(m) {
                                                           colSums(m)
                                                         }))
  summed.counts = summed.counts[order(match(rownames(summed.counts), 
                                            as.character(countsData[, 1]))), ]
  iso.data = getvar(countsData[3:ncol(countsData)], design, 
                    lib.size)
  gene.data = getvar(summed.counts, design, lib.size)
  isoform.numbers = table(table(countsData[, 1]))
  modelString = paste0(functions.string,"data {\n  int<lower=0> Ngenes;\n  int<lower=0> Ntranscripts;\n  int<lower=0> Ncondition1;\n  int<lower=0> Ncondition2;\n  int<lower=2> Nisoforms[Ngenes];\n  real mean_controls[Ngenes];\n  matrix[Ncondition2,Ntranscripts] sd_iso_cases;\n  matrix[Ncondition1,Ntranscripts] sd_iso_controls;\n  matrix[Ncondition2,Ntranscripts] counts_cases;\n  matrix[Ncondition1,Ntranscripts] counts_controls;\n  int<lower=0> num_opt_genes;\n  int run_ids[num_opt_genes];\n  int<lower=1> gene_to_ds[Ngenes];\n\n} \n  parameters { ")
  for (i in (1:length(isoform.numbers))) modelString = paste0(modelString, 
                                                              "\n          \n  simplex[", as.integer(names(isoform.numbers)[i]), 
                                                              "] frac_", names(isoform.numbers)[i], "[", isoform.numbers[i], 
                                                              "];\n                      \n  simplex[", as.integer(names(isoform.numbers)[i]), 
                                                              "] alpha_", names(isoform.numbers)[i], "[", isoform.numbers[i], 
                                                              "];                   \n                      \n                      \n                      ")
  modelString = paste(modelString, "\n  \n  vector[Ngenes] beta;\n  vector[Ngenes] intercept;\n  real<lower=0,upper=1> theta_a;\n  real<lower=0,upper=1> theta_b;\n  \n  }\n  model {\n\n  int l;\n  \n  int idx;  \n\n  int first_i=1;\n\n  int last_i=0;\n\n  ")
  modelString = paste(modelString, "\n  \n  for (genes_itr in run_ids[1:num_opt_genes])\n  {\n    last_i=first_i+Nisoforms[genes_itr]-1;\n\n    intercept[genes_itr] ~ normal(mean_controls[genes_itr],5);\n\n    target+=normal_lpdf_m(beta[genes_itr],theta_b);\n\n  ")
  for (i in names(isoform.numbers)) modelString = paste0(modelString, 
                                                         "if (Nisoforms[genes_itr]==", i, ") {", write.block(i), 
                                                         "};\n ")
  modelString = paste(modelString, "\nfirst_i=last_i+1;\n}\n}")
  stan.mod = stan_model(model_code = modelString)
  tab = iso.data[[1]]
  summed.counts = gene.data[[1]]
  isoform.count = table(countsData[, 1])
  isoform.count = isoform.count[match(unique(countsData[, 
                                                        1]), names(isoform.count))]
  isoform.count = as.numeric(isoform.count)
  split.factor = factor(as.character(countsData[, 1]), levels = unique(as.character(countsData[, 
                                                                                               1])))
  init_l = list(theta_a = 0.5, theta_b = 0.5, beta = rep(0, 
                                                         nrow(summed.counts)), intercept = log2(rowMeans(2^summed.counts[, 
                                                                                                                         labels == 1] - 0.5) + 0.5))
  for (i in as.integer(names(isoform.numbers))) {
    init_l[[length(init_l) + 1]] = matrix(rep(1/i, isoform.numbers[as.integer(names(isoform.numbers)) == 
                                                                     i] * i), ncol = i)
    names(init_l)[length(init_l)] = paste0("alpha_", i)
    row.idxs = countsData[, 1] %in% rownames(summed.counts)[isoform.count == 
                                                              i]
    init_l[[length(init_l) + 1]] = do.call(rbind, lapply(lapply(split(2^tab[row.idxs, 
                                                                            ], countsData[row.idxs, 1]), matrix, ncol = ncol(tab)), 
                                                         function(m) rowSums(m)/sum(rowSums(m))))
    names(init_l)[length(init_l)] = paste0("frac_", i)
  }
  gene.names = unique(countsData[, 1])
  g.group = 1:nrow(summed.counts)
  ds.iso.counts = rep(1, max(isoform.numbers))
  genes.to.ds = c()
  for (i in 1:nrow(summed.counts)) {
    genes.to.ds = c(genes.to.ds, ds.iso.counts[isoform.count[i]])
    ds.iso.counts[isoform.count[i]] = ds.iso.counts[isoform.count[i]] + 
      1
  }
  dataList = list(Ncondition1 = sum(labels == 1), Ncondition2 = sum(labels == 
                                                                      2), counts_cases = t(tab[, labels == 2]), counts_controls = t(tab[, 
                                                                                                                                        labels == 1]), mean_controls = log2(rowMeans(2^summed.counts[, 
                                                                                                                                                                                                     labels == 1] - 0.5) + 0.5), sd_iso_cases = t(sqrt(iso.data[[2]][, 
                                                                                                                                                                                                                                                                     labels == 2])), sd_iso_controls = t(sqrt(iso.data[[2]][, 
                                                                                                                                                                                                                                                                                                                            labels == 1])), Ngenes = nrow(summed.counts), Ntranscripts = nrow(countsData), 
                  Nisoforms = isoform.count, run_ids = g.group, num_opt_genes = length(g.group), 
                  gene_to_ds = genes.to.ds)
  stanFit = optimizing(object = stan.mod, data = dataList, 
                       init = init_l, iter = opt.iter)
  theta_a = stanFit$par[grepl("theta_a", names(stanFit$par))]
  theta_b = stanFit$par[grepl("theta_b", names(stanFit$par))]
  return(list(theta_a, theta_b))
}


write.block2=function (n.iso) 
{
  paste0("\n  idx= gene_to_ds[genes_itr];\n\n  target+=dirichlet_lpdf(frac_", 
         n.iso, "[idx]|rep_vector(1,Nisoforms[genes_itr]));\n  \n  target+=dirichlet_lpdf_m(alpha_", 
         n.iso, "[idx],theta_a);\n  \n  expression_controls[genes_itr,] ~ normal(intercept[genes_itr],sd_exp_controls[genes_itr,]);\n  \n  expression_cases[genes_itr,] ~ normal(intercept[genes_itr]+beta[genes_itr],sd_exp_cases[genes_itr,]);\n\n  l=1;\n  \n  for (j in first_i:last_i)\n  {\n      target+=normal_lpdf(counts_controls[,j] | log2(frac_", 
         n.iso, "[idx][l])+expression_controls[genes_itr,],sd_iso_controls[,j]);\n  \n      target+=normal_lpdf(counts_cases[,j] | log2((frac_", 
         n.iso, "[idx][l]*alpha_", n.iso, "[idx][l])/dot_product(frac_", 
         n.iso, "[idx],alpha_", n.iso, "[idx]))+expression_cases[genes_itr,],sd_iso_cases[,j]);\n    \n      l=l+1;\n  }\n  \n")
}


theta.heirarchical=function (countsData, labels, opt.iter = 30, lib.size = NULL, 
                             n.cores = getOption("mc.cores", 2L)) 
{
  if (sum(labels != labels[order(labels)]) > 0) {
    print("Error: The samples are not ordered!")
    return(NULL)
  }
  labels = factor(labels)
  
  functions.string <- "
  
  functions{
  
  
  real normal_lpdf_m(real y,real theta){
  
  real s=normal_lpdf(y|0,0.04);
  
  real l=normal_lpdf(y|0,5);
  
  if (l>s)
  
  return log(theta)+l;
  
  return log(1-theta)+s;
  
  }
  
  
  real dirichlet_lpdf_m(vector y,real theta){
  
  real s=dirichlet_lpdf(y|rep_vector(50.0,num_elements(y)));
  
  real l=dirichlet_lpdf(y|rep_vector(1.0,num_elements(y)));
  
  if (l>s)
  
  return log(theta)+l;
  
  return log(1-theta)+s;
  
  }
  
  
  
  }\n"

  design = model.matrix(~group, data.frame(group = labels))
  summed.counts = do.call(rbind, lapply(split(countsData[, 
                                                         3:ncol(countsData)], countsData[, 1]), function(m) {
                                                           colSums(m)
                                                         }))
  summed.counts = summed.counts[order(match(rownames(summed.counts), 
                                            as.character(countsData[, 1]))), ]
  iso.data = getvar(countsData[3:ncol(countsData)], design, 
                    lib.size)
  gene.data = getvar(summed.counts, design, lib.size)
  isoform.numbers = table(table(countsData[, 1]))
  modelString = paste0(functions.string,"data {\n  int<lower=0> Ngenes;\n  int<lower=0> Ntranscripts;\n  int<lower=0> Ncondition1;\n  int<lower=0> Ncondition2;\n  int<lower=0> Nisoforms[Ngenes];\n  real mean_controls[Ngenes];\n  matrix[Ngenes,Ncondition2] sd_exp_cases;\n  matrix[Ngenes,Ncondition1] sd_exp_controls;\n  matrix[Ncondition2,Ntranscripts] sd_iso_cases;\n  matrix[Ncondition1,Ntranscripts] sd_iso_controls;\n  matrix[Ncondition2,Ntranscripts] counts_cases;\n  matrix[Ncondition1,Ntranscripts] counts_controls;\n  int<lower=0> num_opt_genes;\n  int run_ids[num_opt_genes];\n  int<lower=1> gene_to_ds[Ngenes];\n  \n}\nparameters {  ")
  for (i in (1:length(isoform.numbers))) modelString = paste0(modelString, 
                                                              "\n                       \n        simplex[", as.integer(names(isoform.numbers)[i]), 
                                                              "] frac_", names(isoform.numbers)[i], "[", isoform.numbers[i], 
                                                              "];\n                       \n        simplex[", as.integer(names(isoform.numbers)[i]), 
                                                              "] alpha_", names(isoform.numbers)[i], "[", isoform.numbers[i], 
                                                              "];                   \n                       \n            ")
  modelString = paste(modelString, "\n\n  vector[Ngenes] beta;\n  vector[Ngenes] intercept;\n  vector<lower=0,upper=1>[Ntranscripts] alpha;\n  vector<lower=0,upper=1>[Ntranscripts] frac;\n  matrix[Ngenes,Ncondition1] expression_controls;\n  matrix[Ngenes,Ncondition2] expression_cases;\n  real<lower=0,upper=1> theta_a;\n  real<lower=0,upper=1> theta_b;\n    \n  }\n  model {\n \n  int l;\n  \n  int idx;  \n\n  int first_i=1;\n\n  int last_i=0;\n  \n  for (genes_itr in run_ids[1:num_opt_genes])\n  {\n      last_i=first_i+Nisoforms[genes_itr]-1;\n  \n      target+=normal_lpdf_m(beta[genes_itr],theta_b);\n  \n      intercept[genes_itr] ~ normal(mean_controls[genes_itr],5);\n ")
  for (i in names(isoform.numbers)) modelString = paste0(modelString, 
                                                         "if (Nisoforms[genes_itr]==", i, ") {", write.block2(i), 
                                                         "};\n ")
  modelString = paste(modelString, "\nfirst_i=last_i+1;\n}\n}")
  stan.mod = stan_model(model_code = modelString)
  tab = iso.data[[1]]
  summed.counts = gene.data[[1]]
  isoform.count = table(countsData[, 1])
  isoform.count = isoform.count[match(unique(countsData[, 
                                                        1]), names(isoform.count))]
  isoform.count = as.numeric(isoform.count)
  split.factor = factor(as.character(countsData[, 1]), levels = unique(as.character(countsData[, 
                                                                                               1])))
  init_l = list(theta_a = 0.5, theta_b = 0.5, beta = rep(0, 
                                                         nrow(summed.counts)), intercept = log2(rowMeans(2^summed.counts[, 
                                                                                                                         labels == 1] - 0.5) + 0.5), expression_cases = summed.counts[, 
                                                                                                                                                                                      labels == 2], expression_controls = summed.counts[, 
                                                                                                                                                                                                                                        labels == 1])
  for (i in as.integer(names(isoform.numbers))) {
    init_l[[length(init_l) + 1]] = matrix(rep(1/i, isoform.numbers[as.integer(names(isoform.numbers)) == 
                                                                     i] * i), ncol = i)
    names(init_l)[length(init_l)] = paste0("alpha_", i)
    row.idxs = countsData[, 1] %in% rownames(summed.counts)[isoform.count == 
                                                              i]
    init_l[[length(init_l) + 1]] = do.call(rbind, lapply(lapply(split(2^tab[row.idxs, 
                                                                            ], countsData[row.idxs, 1]), matrix, ncol = ncol(tab)), 
                                                         function(m) rowSums(m)/sum(rowSums(m))))
    names(init_l)[length(init_l)] = paste0("frac_", i)
  }
  gene.names = unique(countsData[, 1])
  g.group = 1:nrow(summed.counts)
  ds.iso.counts = rep(1, max(isoform.numbers))
  genes.to.ds = c()
  for (i in 1:nrow(summed.counts)) {
    genes.to.ds = c(genes.to.ds, ds.iso.counts[isoform.count[i]])
    ds.iso.counts[isoform.count[i]] = ds.iso.counts[isoform.count[i]] + 
      1
  }
  dataList = list(Ncondition1 = sum(labels == 1), Ncondition2 = sum(labels == 
                                                                      2), counts_cases = t(tab[, labels == 2]), counts_controls = t(tab[, 
                                                                                                                                        labels == 1]), mean_controls = log2(rowMeans(2^summed.counts[, 
                                                                                                                                                                                                     labels == 1] - 0.5) + 0.5), sd_iso_cases = t(sqrt(iso.data[[2]][, 
                                                                                                                                                                                                                                                                     labels == 2])), sd_iso_controls = t(sqrt(iso.data[[2]][, 
                                                                                                                                                                                                                                                                                                                            labels == 1])), sd_exp_cases = sqrt(gene.data[[2]][, 
                                                                                                                                                                                                                                                                                                                                                                               labels == 2]), sd_exp_controls = sqrt(gene.data[[2]][, 
                                                                                                                                                                                                                                                                                                                                                                                                                                    labels == 1]), Ngenes = nrow(summed.counts), Ntranscripts = nrow(countsData), 
                  Nisoforms = isoform.count, run_ids = g.group, num_opt_genes = length(g.group), 
                  gene_to_ds = genes.to.ds)
  stanFit = optimizing(object = stan.mod, data = dataList, 
                       init = init_l, iter = opt.iter)
  theta_a = stanFit$par[grepl("theta_a", names(stanFit$par))]
  theta_b = stanFit$par[grepl("theta_b", names(stanFit$par))]
  return(list(theta_a, theta_b))
}
