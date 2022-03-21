.. _rsttranscriptome:

===============================
Inferring transcript readcounts
===============================

For this tutorial, we will use the  isoform quantification tool
`RSEM <https://pubmed.ncbi.nlm.nih.gov/21816040/>`_ to generate isoform counts.
A typical run of RSEM consists of just two steps. First, a set of reference transcript 
sequences are generated and preprocessed for use by later RSEM steps. Second, a set of 
RNA-Seq reads are aligned to the reference transcripts and the resulting alignments are 
used to estimate abundances and their credibility intervals. 



Installing rsem
^^^^^^^^^^^^^^^



Preparing the reference genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, run the following commands to download and unzip the human genome and gene annotation:


.. code-block:: bash

    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.p13.genome.fa.gz
    gunzip GRCh38.p13.genome.fa.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
    gunzip gencode.v39.annotation.gtf.gz


Create a reference using RSEM by executing the following command:

.. code-block:: bash
    
    rsem-prepare-reference --gtf gencode.v39.annotation.gtf --num-threads 8 --star GRCh38.p13.genome.fa ref_rsem/ref_rsem