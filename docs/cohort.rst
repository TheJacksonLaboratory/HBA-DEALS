.. _rstcohort:

===================
Cohort RNA-seq data
===================


HBA-DEALS is designed for the analysis of RNA-seq cohort data. The nature of 
the cohorts is arbitrary but will commonly be labeled as cases and controls or group 1 vs. group2.


In general, there are many ways to obtain such data including generating RNA-seq data or downloading it from an appropriate site. 
For this tutorial, we download a dataset from the `NCBI Sequence Read Archive (SRA) <https://www.ncbi.nlm.nih.gov/sra>`_, 
which is the primary archive of NGS datasets.


Downloading from SRA
^^^^^^^^^^^^^^^^^^^^

We will use the `SRA Toolkit <https://hpc.nih.gov/apps/sratoolkit.html>`_.

We will use 8 RNA-seq samples from the dataset  `SRP149366 <https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP149366>`_,
which investigated estrogen responsive transcriptome of estrogen receptor positive normal human breast cells in 3D cultures.


+-------------+-------------+
| control     | estradiol   |
+=============+=============+
| SRR7236472  |  SRR7236480 |
| SRR7236473  |  SRR7236481 |
| SRR7236474  |  SRR7236482 |
| SRR7236475  |  SRR7236483 |
+-------------+-------------+





First install ``fasterq-dump`` on your system accorrding to the SRA toolkit instructions. The `download page <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`_ 
may be helpful.

Then execute the following command from
the shell.


.. code-block:: bash
   :caption: Downloading RNA-seq files with the SRA Toolkit

    for srr in SRR7236472 SRR7236473 SRR7236474 SRR7236475 SRR7236480 SRR7236481 SRR7236482 SRR7236483; do \
        prefetch $srr 
        fasterq-dump -t tmp/ --split-files --threads 8 --outdir tutorial/ $srr 
    done

If necessary, change the ``--threads`` argument according to the resources of your system. 
The downloaded ``*.fastq`` files will be written to the ``tutorial`` directory. 
Other directories that were created (e.g., ``SRR7236472``, which contains the file ``SRR7236472.sra``) can be deleted.

This step will typically take a few hours.


Cleaning the reads
^^^^^^^^^^^^^^^^^^

`fastp <https://academic.oup.com/bioinformatics/article/34/17/i884/5093234>`_ performs
quality control, adapter trimming, quality filtering, per-read quality pruning and many other 
operations with a single scan of the FASTQ data.

fastp can be downloaded at its `GitHub site <https://github.com/OpenGene/fastp>`_, which also has installation 
instructions for various platforms. If you have a debian or Ubuntu system, then the easiest installation is just

.. code-block:: bash

    sudo apt install fastp



After you have installed fastp, run the following command in the shell.


.. code-block:: bash
   :caption: Cleaning FASTQ files with fastp

    for srr in SRR7236472 SRR7236473 SRR7236474 SRR7236475 SRR7236480 SRR7236481 SRR7236482 SRR7236483; do \
        fastp -i tutorial/${srr}_1.fastq -I tutorial/${srr}_2.fastq -o tutorial/${srr}_trimmed_1.fastq -O tutorial/${srr}_trimmed_2.fastq; \
    done

This will create for each .fastq file that was downloaded a file by the same name with an added "_trimmed_" in its name. 
At this point, you can delete the original ``*.fastq`` files if desired.




