.. _rsthbadeals:

=================
Running HBA-DEALS
=================



We will now read the expected counts from the rsem directory and create a count matrix, which we will pass as input to HBA-DEALS.

The following code creates the counts matrix


.. code-block:: R 

    countsData=NULL  #The counts matrix
    labels=c()   #The sample labels, 1 for control and 2 for case
    for (srr in c('SRR7236472','SRR7236473','SRR7236474','SRR7236475')) {
        if (is.null(countsData))
            countsData=read.table(paste0('rsem/',srr,'.isoforms.results'),header=TRUE)[c(2,1)]
        next.file=read.table(paste0('rsem/',srr,'.isoforms.results'),header=TRUE)
        countsData=cbind(countsData,next.file$expected_count)
       labels=c(labels,1)
    }
    for (srr in c('SRR7236480','SRR7236481','SRR7236482','SRR7236483')) {
       next.file=read.table(paste0('rsem/',srr,'.isoforms.results'),header=TRUE)
       countsData=cbind(countsData,next.file$expected_count)
       labels=c(labels,2)
    }
    countsData=countsData[rowSums(countsData[,-c(1,2)]>0)>=ncol(countsData)-2,] #remove low-count isoforms
    num.iso=unlist(lapply(countsData$gene_id,function(x){sum(countsData$gene_id %in% x)}))
    countsData=countsData[num.iso>1,]  #keep only genes with at least two isoforms


The following code performs the HBA-DEALS analysis.

.. code-block:: R 

    library(hbadeals)
    library(edgeR)
    res=hbadeals(countsData = countsData,labels = labels,n.cores = 32,isoform.level=TRUE,mcmc.iter=100000,mcmc.warmup=10000,mtc=TRUE,
             lib.size=colSums(countsData[,-c(1,2)])*calcNormFactors(as.matrix(countsData[,-c(1,2)]),method='TMM'))
    write.table(res,'hbadeals_output.txt', sep='\t',quote = F,col.names = T,row.names = F)