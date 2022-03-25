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




			
HBA-DEALS arguments
^^^^^^^^^^^^^^^^^^^

* countsData: A table of gene names, transcript names and transcript counts in each sample.  At least two transcripts must correspond to
   each gene.
* labels:  An ordered vector of 1's and 2's.  Its length is ncol(countsdata)-2.  Each entry indicates whether the corresponding sample/column
  of countsData belongs to the first experimental condition or the second.
* n.cores: The number of cores to use in the calculation.  It is recommended to dedicate as many cores as possible.
* isoform.level: if ``true``, return 1-probability of differential proportion for each transcript be returned. If ``FALSE``  (the default) return 1-probability of differential splicing for each gene.
* mcmc.iter: The number of iterations of the MCMC algorithm after warmup.
* mcmc.warmup: The number of warmup iterations of the MCMC algorithm.
* hierarchy:  Determines whether a hierarchical model will be used (hierarchy='yes'), a flat one (hierarchy='no') or will the decision
  will be made automatically (hierarchy='auto')
* lib.size:  A numeric vector containing total library sizes for each sample. If not provided, the default is columnwise count totals.
* mtc:  A logical argument (default FALSE) that indicates whether the output probabilities should be corrected for multiple comparisons.


Output format
^^^^^^^^^^^^^

The HBA-DEALS output file contains 4 columns. The first column is the gene name, the second is the transcript name,
the third is the fold change, and the fourth is 1-probability of differential expression or proportion(splicing),
which is the posterior error probability (PEP).  Entries that
refer to expression have 'Expression' in their second column.  If isoform.level is FALSE, entries that refer to differential 
splicing of the gene will have 'Splicing' in their second column entry.  The fold change for expression is 
given as log2 fold change, and for splicing as fold change.




+-------------------+-------------------+-------------------+---- -----+
| Gene              | Isoform           | ExplogFC/FC       |P         | 
+===================+===================+===================+==========+
| ENSG00000000419	| Expression	    |0.00360856329234098|	0.98668|
+-------------------+-------------------+-------------------+---- -----+
| ENSG00000000419	| ENST00000371588	| 0.996357449743301	| 0.95506  |
+-------------------+-------------------+-------------------+---- -----+
| ENSG00000196141	| Expression	    |-2.74059622402359	| 0        |
+-------------------+-------------------+-------------------+---- -----+
| ENSG00000196233	| ENST00000371103	|0.10984085543658	| 0.04972  |
+-------------------+-------------------+-------------------+---- -----+	

