hbadeals.heirarchy=function (countsData, labels, n.cores = getOption("mc.cores", 
                                                                     2L), isoform.level = FALSE, mcmc.iter = 3000, mcmc.warmup = 4000, 
                             lib.size = NULL, mtc = FALSE) 
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
  if (mtc) {
    theta.vals = theta.heirarchical(countsData, labels, 
                                    1000, lib.size, n.cores)
    theta_a = theta.vals[[1]]
    theta_b = theta.vals[[2]]
  }
  modelString = paste0(functions.string,"\n  data {\n    int<lower=0> Nisoforms;\n    int<lower=0> Ncondition1;\n    int<lower=0> Ncondition2;\n    real mean_controls;\n    vector[Ncondition1] sd_exp_controls;\n    vector[Ncondition2] sd_exp_cases;\n    vector[Nisoforms] sd_iso_cases[Ncondition2];\n    vector[Nisoforms] sd_iso_controls[Ncondition1];\n    vector[Nisoforms] counts_cases[Ncondition2];\n    vector[Nisoforms] counts_controls[Ncondition1];\n    }\n\n  parameters {\n    simplex[Nisoforms] frac;\n    vector[Ncondition1] expression_controls;\n    vector[Ncondition2] expression_cases;\n    real beta;\n    real intercept;\n    simplex[Nisoforms] alpha;\n    }\n\n  model {")
  if (mtc) {
    modelString = paste0(modelString, "target+=normal_lpdf_m(beta,",theta_b, 
                         ");\n\n                       \n                                     target+=dirichlet_lpdf_m(alpha,", 
                         theta_a, ");")
  }
  else {
    modelString = paste0(modelString, "alpha ~ dirichlet(rep_vector(1.0,Nisoforms));\nbeta ~ normal(0,5);\n")
  }
  modelString = paste0(modelString, "\n\n      frac ~ dirichlet(rep_vector(1.0,Nisoforms));\n\n      intercept ~ normal(mean_controls,5);\n\n      expression_controls ~ normal(intercept,sd_exp_controls);\n\n      expression_cases ~ normal(intercept+beta,sd_exp_cases);\n\n      for ( j in 1:Nisoforms )\n      {\n          counts_controls[,j] ~ normal(log2(frac[j])+expression_controls,sd_iso_controls[,j]);\n\n          counts_cases[,j] ~ normal(log2((frac[j]*alpha[j])/dot_product(frac,alpha))+expression_cases,sd_iso_cases[,j]);\n      }\n\n  }\n  ")
  stan.mod = stan_model(model_code = modelString)
  tab = iso.data[[1]]
  summed.counts = gene.data[[1]]
  results.tab = do.call(rbind, mclapply(1:length(unique(countsData[, 
                                                                   1])), function(gene.number) {
                                                                     gene.rows = which(countsData[, 1] == unique(countsData[, 
                                                                                                                            1])[gene.number])
                                                                     num.isoforms = length(gene.rows)
                                                                     initf = function() {
                                                                       list(frac = rowSums(2^tab[gene.rows, ])/sum(rowSums(2^tab[gene.rows, 
                                                                                                                                 ])), expression_cases = summed.counts[gene.number, 
                                                                                                                                                                       labels == 2], expression_controls = summed.counts[gene.number, 
                                                                                                                                                                                                                         labels == 1], alpha = array(rep(1/num.isoforms, 
                                                                                                                                                                                                                                                         num.isoforms), dim = num.isoforms), beta = 0, 
                                                                            intercept = log2(mean(2^summed.counts[gene.number, 
                                                                                                                  labels == 1] - 0.5) + 0.5))
                                                                     }
                                                                     res = matrix(ncol = 4, nrow = 0)
                                                                     colnames(res) = c("Gene", "Isoform", "ExplogFC/FC", 
                                                                                       "P")
                                                                     dataList = list(counts_cases = t(tab[gene.rows, labels == 
                                                                                                            2]), counts_controls = t(tab[gene.rows, labels == 
                                                                                                                                           1]), Nisoforms = num.isoforms, Ncondition1 = sum(labels == 
                                                                                                                                                                                              1), Ncondition2 = sum(labels == 2), mean_controls = log2(mean(2^summed.counts[gene.number, 
                                                                                                                                                                                                                                                                            labels == 1] - 0.5) + 0.5), sd_iso_cases = t(sqrt(iso.data[[2]][gene.rows, 
                                                                                                                                                                                                                                                                                                                                            labels == 2])), sd_iso_controls = t(sqrt(iso.data[[2]][gene.rows, 
                                                                                                                                                                                                                                                                                                                                                                                                   labels == 1])), sd_exp_cases = sqrt(gene.data[[2]][gene.number, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                      labels == 2]), sd_exp_controls = sqrt(gene.data[[2]][gene.number, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           labels == 1]))
                                                                     stanFit = sampling(object = stan.mod, data = dataList, 
                                                                                        cores = 1, init = initf, chains = 1, refresh = 0, 
                                                                                        iter = mcmc.warmup + mcmc.iter, warmup = mcmc.warmup, 
                                                                                        thin = 1)
                                                                     mcmcCoda = mcmc.list(lapply(1:ncol(stanFit), function(x) {
                                                                       mcmc(as.array(stanFit)[, x, ])
                                                                     }))
                                                                     beta.col = which(colnames(mcmcCoda[[1]]) == "beta")
                                                                     mean.de = mean(mcmcCoda[[1]][, beta.col])
                                                                     expression.p = p_rope(as.numeric(mcmcCoda[[1]][, beta.col]), 
                                                                                           range = c(-0.1, 0.1))$p_ROPE
                                                                     gene.name = as.character(unique(countsData[, 1])[gene.number])
                                                                     res = rbind(res, c(gene.name, "Expression", mean.de, 
                                                                                        expression.p))
                                                                     frac = rep(0, num.isoforms)
                                                                     alpha = rep(0, num.isoforms)
                                                                     alpha.p = rep(0, num.isoforms)
                                                                     for (i in (1:num.isoforms)) {
                                                                       next.frac.col = which(colnames(mcmcCoda[[1]]) == 
                                                                                               paste0("frac[", i, "]"))
                                                                       frac[i] = mean(mcmcCoda[[1]][, next.frac.col])
                                                                       next.alpha.col = which(colnames(mcmcCoda[[1]]) == 
                                                                                                paste0("alpha[", i, "]"))
                                                                       alpha[i] = mean(mcmcCoda[[1]][, next.alpha.col])
                                                                     }
                                                                     if (isoform.level) {
                                                                       uni.val = sum(frac * alpha)
                                                                     }
                                                                     else {
                                                                       uni.val = 1/num.isoforms
                                                                     }
                                                                     for (i in (1:num.isoforms)) {
                                                                       next.alpha.col = which(colnames(mcmcCoda[[1]]) == 
                                                                                                paste0("alpha[", i, "]"))
                                                                       alpha.p[i] = p_rope(as.numeric(mcmcCoda[[1]][, next.alpha.col]) - 
                                                                                             uni.val, range = c(-0.2 * uni.val, 0.2 * uni.val))$p_ROPE
                                                                     }
                                                                     frac.2 = (frac * alpha)/sum(frac * alpha)
                                                                     fc.ds = frac.2/frac
                                                                     if (isoform.level) {
                                                                       for (i in (1:num.isoforms)) res = rbind(res, c(gene.name, 
                                                                                                                      as.character(countsData[gene.rows[i], 2]), fc.ds[i], 
                                                                                                                      alpha.p[i]))
                                                                     }
                                                                     else {
                                                                       min.p = min(alpha.p)
                                                                       res = rbind(res, c(gene.name, "Splicing", fc.ds[which(alpha.p == 
                                                                                                                               min(alpha.p))[1]], min.p))
                                                                     }
                                                                     return(res)
                                                                   }, mc.cores = n.cores))
  return(results.tab)
}



hbadeals.flat=function (countsData, labels, n.cores = getOption("mc.cores", 
                                                                2L), isoform.level = FALSE, mcmc.iter = 3000, mcmc.warmup = 4000, 
                        lib.size = NULL, mtc = TRUE) 
{
  if (sum(labels != labels[order(labels)]) > 0) {
    print("Error: The samples are not ordered!")
    return(NULL)
  }
  
  
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
  
  
  labels = factor(labels)
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
  if (mtc) {
    theta.vals = theta.flat(countsData, labels, 1000, lib.size, 
                            n.cores)
    theta_a = theta.vals[[1]]
    theta_b = theta.vals[[2]]
  }
  modelString = paste0(functions.string,"\n    data {\n        int<lower=0> Nisoforms;\n        int<lower=0> Ncondition1;\n        int<lower=0> Ncondition2;\n        real mean_controls;\n        vector[Nisoforms] sd_iso_cases[Ncondition2];\n        vector[Nisoforms] sd_iso_controls[Ncondition1];\n        vector[Nisoforms] counts_cases[Ncondition2];\n        vector[Nisoforms] counts_controls[Ncondition1];\n    }\n    \n    parameters {\n        simplex[Nisoforms] frac;\n        real beta;\n        real intercept;\n        simplex[Nisoforms] alpha;\n    }\n    \n    model {")
  if (mtc) {
    modelString = paste0(modelString, "target+=normal_lpdf_m(beta,",theta_b, 
                         ")+dirichlet_lpdf_m(alpha,", 
                         theta_a, ");")
  }
  else {
    modelString = paste0(modelString, "alpha ~ dirichlet(rep_vector(1.0,Nisoforms));\nbeta ~ normal(0,5);\n")
  }
  modelString = paste0(modelString, "\n        \n        frac ~ dirichlet(rep_vector(1.0,Nisoforms));\n      \n        intercept ~ normal(mean_controls,5);\n        \n        for ( j in 1:Nisoforms )\n        {\n            counts_controls[,j] ~ normal(log2(frac[j])+intercept,sd_iso_controls[,j]);\n            \n            counts_cases[,j] ~ normal(log2((frac[j]*alpha[j])/dot_product(frac,alpha))+intercept+beta,sd_iso_cases[,j]);\n        }\n        \n    }\n    ")
  stan.mod = stan_model(model_code = modelString)
  tab = iso.data[[1]]
  summed.counts = gene.data[[1]]
  results.tab = do.call(rbind, mclapply(1:length(unique(countsData[, 
                                                                   1])), function(gene.number) {
                                                                     gene.rows = which(countsData[, 1] == unique(countsData[, 
                                                                                                                            1])[gene.number])
                                                                     num.isoforms = length(gene.rows)
                                                                     initf = function() {
                                                                       list(frac = rowSums(2^tab[gene.rows, ])/sum(rowSums(2^tab[gene.rows, 
                                                                                                                                 ])), alpha = array(rep(1/num.isoforms, num.isoforms), 
                                                                                                                                                    dim = num.isoforms), beta = 0, intercept = log2(mean(2^summed.counts[gene.number, 
                                                                                                                                                                                                                         labels == 1] - 0.5) + 0.5))
                                                                     }
                                                                     res = matrix(ncol = 4, nrow = 0)
                                                                     colnames(res) = c("Gene", "Isoform", "ExplogFC/FC", 
                                                                                       "P")
                                                                     dataList = list(counts_cases = t(tab[gene.rows, labels == 
                                                                                                            2]), counts_controls = t(tab[gene.rows, labels == 
                                                                                                                                           1]), Nisoforms = num.isoforms, Ncondition1 = sum(labels == 
                                                                                                                                                                                              1), Ncondition2 = sum(labels == 2), mean_controls = log2(mean(2^summed.counts[gene.number, 
                                                                                                                                                                                                                                                                            labels == 1] - 0.5) + 0.5), sd_iso_cases = t(sqrt(iso.data[[2]][gene.rows, 
                                                                                                                                                                                                                                                                                                                                            labels == 2])), sd_iso_controls = t(sqrt(iso.data[[2]][gene.rows, 
                                                                                                                                                                                                                                                                                                                                                                                                   labels == 1])))
                                                                     stanFit = sampling(object = stan.mod, data = dataList, 
                                                                                        cores = 1, init = initf, chains = 1, refresh = 0, 
                                                                                        iter = mcmc.warmup + mcmc.iter, warmup = mcmc.warmup, 
                                                                                        thin = 1)
                                                                     mcmcCoda = mcmc.list(lapply(1:ncol(stanFit), function(x) {
                                                                       mcmc(as.array(stanFit)[, x, ])
                                                                     }))
                                                                     beta.col = which(colnames(mcmcCoda[[1]]) == "beta")
                                                                     mean.de = mean(mcmcCoda[[1]][, beta.col])
                                                                     expression.p = p_rope(as.numeric(mcmcCoda[[1]][, beta.col]), 
                                                                                           range = c(-0.1, 0.1))$p_ROPE
                                                                     gene.name = as.character(unique(countsData[, 1])[gene.number])
                                                                     res = rbind(res, c(gene.name, "Expression", mean.de, 
                                                                                        expression.p))
                                                                     frac = rep(0, num.isoforms)
                                                                     alpha = rep(0, num.isoforms)
                                                                     alpha.p = rep(0, num.isoforms)
                                                                     for (i in (1:num.isoforms)) {
                                                                       next.frac.col = which(colnames(mcmcCoda[[1]]) == 
                                                                                               paste0("frac[", i, "]"))
                                                                       frac[i] = mean(mcmcCoda[[1]][, next.frac.col])
                                                                       next.alpha.col = which(colnames(mcmcCoda[[1]]) == 
                                                                                                paste0("alpha[", i, "]"))
                                                                       alpha[i] = mean(mcmcCoda[[1]][, next.alpha.col])
                                                                     }
                                                                     if (isoform.level) {
                                                                       uni.val = sum(frac * alpha)
                                                                     }
                                                                     else {
                                                                       uni.val = 1/num.isoforms
                                                                     }
                                                                     for (i in (1:num.isoforms)) {
                                                                       next.alpha.col = which(colnames(mcmcCoda[[1]]) == 
                                                                                                paste0("alpha[", i, "]"))
                                                                       alpha.p[i] = p_rope(as.numeric(mcmcCoda[[1]][, next.alpha.col]) - 
                                                                                             uni.val, range = c(-0.2 * uni.val, 0.2 * uni.val))$p_ROPE
                                                                     }
                                                                     frac.2 = (frac * alpha)/sum(frac * alpha)
                                                                     fc.ds = frac.2/frac
                                                                     if (isoform.level) {
                                                                       for (i in (1:num.isoforms)) res = rbind(res, c(gene.name, 
                                                                                                                      as.character(countsData[gene.rows[i], 2]), fc.ds[i], 
                                                                                                                      alpha.p[i]))
                                                                     }
                                                                     else {
                                                                       min.p = min(alpha.p)
                                                                       res = rbind(res, c(gene.name, "Splicing", fc.ds[which(alpha.p == 
                                                                                                                               min(alpha.p))[1]], min.p))
                                                                     }
                                                                     return(res)
                                                                   }, mc.cores = n.cores))
  return(results.tab)
}


hbadeals=function (countsData, labels, n.cores = getOption("mc.cores", 
                                                           2L), isoform.level = FALSE, mcmc.iter = 3000, mcmc.warmup = 4000, 
                   hierarchy = "no", lib.size = NULL, mtc = TRUE) 
{
  use.heirarchical = TRUE
  if (hierarchy == "auto") {
    set.seed(123)
    use.heirarchical = choose.model(countsData = countsData, 
                                    labels = labels, n.cores = n.cores, num.genes = 500, 
                                    mcmc.iter = mcmc.iter, mcmc.warmup = mcmc.warmup, 
                                    lib.size = lib.size)
  }
  else if (hierarchy == "no") {
    use.heirarchical = FALSE
  }
  if (use.heirarchical) 
    return(hbadeals.heirarchy(countsData = countsData, labels = labels, 
                              n.cores = n.cores, isoform.level = isoform.level, 
                              mcmc.iter = mcmc.iter, mcmc.warmup = mcmc.warmup, 
                              lib.size = lib.size, mtc = mtc))
  return(hbadeals.flat(countsData = countsData, labels = labels, 
                       n.cores = n.cores, isoform.level = isoform.level, mcmc.iter = mcmc.iter, 
                       mcmc.warmup = mcmc.warmup, lib.size = lib.size, mtc = mtc))
}
