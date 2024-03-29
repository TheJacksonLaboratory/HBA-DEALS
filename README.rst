##############################################################################################
Hierarchical Bayesian Analysis of Differential Expression and ALternative Splicing (HBA-DEALS)
##############################################################################################

HBA-DEALS is an R package that simultaneously characterizes differential gene expression and 
alternative splicing in high-throughput gene expression data.  It uses counts data for isoforms (alternative transcripts of a gene)
to infer the parameters of a hierarchical Bayesian model of expression and splicing.  It then uses the posterior of the parameters to determine the existence of differential expression and/or differential alternative splicing.  Isoform counts can be derived from short-read sequencing RNA-Seq data or from long-read RNA sequencing data such as that generated by the PacBio SequelII platform.

A manuscript describing HBA-DEALS provides details on the algorithm: `Karlebach et al, 2020, Genome Biology 21:171 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02072-6>`_



In order to install the package, in R use the command:  ::

  devtools::install_github("TheJacksonLaboratory/HBA-DEALS")

And then load the library with: ::

  library(hbadeals)

The package provides two functions.  The function 'simulate' that generates a sample dataset of isoform counts.  The function's output appears in the working directory.  The file 'counts.txt' contains counts in the input format that the tool expects.  The files 'de_genes.txt' and 'ds_isoforms.txt' contain the identifiers of genes and isoforms that were differentially expressed and spliced, respectively.
The second function, 'hbadeals', takes as input a data frame of counts and a vector that contains labels for cases and control samples, and computes the probabilities of differential expression and splicing.  In order to use this function, the counts data has to be read into a data frame and passed to the function as an argument.  It is recommended to use the parameter mcmc.cores in order to allocate as many cores as possible to the execution of the function, as the MCMC algorithm that HBA-DEALS is computationally demanding.  Increasing the value of the parameter mcmc.iter will increase the accuracy of the probabilities computed by HBA-DEALS, as it corresponds to the number of MCMC steps.  Similarly, mcmc.warmup is the number of warmup iterations and can affect accuracy if too small.  The most effective way to use HBA-DEALS is as part of an NGS pipeline rather than as a desktop application.

The package has dependencies that should be installed with it.  Please see the documentation of the relevant packages if encountering any platform-specific issues.
