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

Download the ``RSEM-1.3.3.tar.gz`` file, uncompress it, and build it.

.. code-block:: bash

    tar xfvz RSEM-1.3.3.tar.gz
    cd RSEM-1.3.3
    make

You can follow this with ``sudo make install``, softlink to the required executables, or add the path of this directory to the path
when you execute the RSEM commands below. The latter option is as follows.

.. code-block:: bash

    export PATH=/home/<user>/<somepath>/RSEM-1.3.3/:$PATH


Installing STAR
^^^^^^^^^^^^^^^

RSEM requires `Spliced Transcripts Alignment to a Reference (STAR) <https://pubmed.ncbi.nlm.nih.gov/23104886/>`_
to be installed and available in the path. STAR can be downloaded from its `GitHub site <https://github.com/alexdobin/STAR>`_. 
For instance, download the statically compiled executable (``STAR_Linux_x86_64_static.zip``), unpack it, and put it in the PATH similarly to the above.





Preparing the reference genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, run the following commands to download and unzip the human genome and gene annotation:


.. code-block:: bash

    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.p13.genome.fa.gz
    gunzip GRCh38.p13.genome.fa.gz
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
    gunzip gencode.v39.annotation.gtf.gz


Create a reference using RSEM by executing the following commands (note that we make a directory for the output of ths RSEM script).

.. code-block:: bash
    
    mkdir ref_rsem
    rsem-prepare-reference --gtf gencode.v39.annotation.gtf --num-threads 8 --star GRCh38.p13.genome.fa ref_rsem/ref_rsem


Note that the RSEM steps are computationally intensive and would typically be done in a high-performance computing environment rather than 
on a laptop.

Counting isoforms
^^^^^^^^^^^^^^^^^

RSEM is a software package for estimating gene and isoform expression levels from RNA-Seq data. 

.. code-block:: bash

    for srr in SRR7236472 SRR7236473 SRR7236474 SRR7236475 SRR7236480 SRR7236481 SRR7236482 SRR7236483; do 
        rsem-calculate-expression --star --paired-end --no-bam-output --num-threads 8 tutorial/${srr}_trimmed_1.fastq tutorial/${srr}_trimmed_2.fastq ref_rsem/ref_rsem rsem/$srr; 
    done



This will produce outfiles with names such as ``SRR7236472.isoforms.results`` with the following table structure.
The HBA-DEALS data ingest script will use this file as input.


+-------------------+-------------------+-------------+-----------------+----------------+---------+---------+---------+
| transcript_id     | gene_id           | length      |effective_length | expected_count |  TPM    |  FPKM   | IsoPct  |
+===================+===================+=============+=================+================+=========+=========+=========+
| ENST00000373020.9 | ENSG00000000003.15|3768         | 3546.60         |  1265.04       |  26.44  | 20.76   | 77.37   +
+-------------------+-------------------+-------------+-----------------+----------------+---------+---------+---------+
| ENST00000494424.1 | ENSG00000000003.15|820          | 598.61          |  0.00          |  0.00   |  0.00   | 0.00    +
+-------------------+-------------------+-------------+-----------------+----------------+---------+---------+---------+
| ENST00000496771.5 | ENSG00000000003.15|1025         | 803.60          |  20.42         |  1.88   |  1.48   | 5.51    +
+-------------------+-------------------+-------------+-----------------+----------------+---------+---------+---------+



                


