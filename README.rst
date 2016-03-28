=======
ezBAMQC
=======

*"ezBAMQC, a tool to check the quality of mapped next generation sequencing files."*

:Codeology Icon:

   .. image:: https://raw.githubusercontent.com/mhammell-laboratory/bamqc/master/doc/bamqc-icon.gif
     :alt: generated at codeology.braintreepayments.com/mhammell-laboratory/bamqc
     :target: http://codeology.braintreepayments.com/mhammell-laboratory/bamqc

:Description:

   ezBAMQC is a tool to check the quality of either one or many mapped next-generation-sequencing datasets. It conducts comprehensive evaluations of aligned sequencing data from multiple aspects including: clipping profile, mapping quality distribution, mapped read length distribution, genomic/transcriptomic mapping distribution, inner distance distribution (for paired-end reads), ribosomal RNA contamination, transcript 5’ and 3’ end bias, transcription dropout rate, sample correlations, sample reproducibility, sample variations. It outputs a set of tables and plots and one HTML page that contains a summary of the results. Many metrics are designed for RNA-seq data specifically, but ezBAMQC can be applied to any mapped sequencing dataset such as RNA-seq, CLIP-seq, GRO-seq, ChIP-seq, DNA-seq and so on.

:Links:

    `Github Page <https://github.com/mhammell-laboratory/bamqc>`_

    `Pypi Page <https://pypi.python.org/pypi/ezBAMQC>`_

    `MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

:Authors:
    Ying Jin, David Molik, and Molly Hammell

:Version: 0.6.7

:Contact:
    Ying Jin (yjin@cshl.edu)

Installation guide for ezBAMQC for from source installs
=======================================================

When installing ezBAMQC there are several options, but the main point is: since ezBAMQC uses C++ STD 11 you'll need a version of GCC that can support that, this useally means 4.8 or 4.9. beyond that, you'll need Python, R and Corrplot for interfacing with the C code.

:Intallation:
   `Source Code <https://github.com/mhammell-laboratory/ezBAMQC/releases>`_

   `Pypi <https://pypi.python.org/pypi?:action=display&name=ezBAMQC>`_

:Prerequisites:
    * `python2.7 <https://www.python.org/download/releases/2.7/>`_
    * `R <https://www.r-project.org/>`_
    * `corrplot <https://cran.r-project.org/web/packages/corrplot/>`_
    * `GCC 4.8.1 or greater <https://gcc.gnu.org/gcc-4.8/>`_ GCC 4.9.1 or greater is recomended for PyPi install 

:Notes:
    * While there are multiple methods of installing the prerequistes it may help to look at (if using a yum based linux distro):*
    * `Devtoolset-3 <https://access.redhat.com/documentation/en-US/Red_Hat_Developer_Toolset/3/html/User_Guide/sect-Red_Hat_Developer_Toolset-Install.html>`_ for GCC compilers
    * `IUS <https://ius.io/>`_ for Python2.7
    * `Software Collections <https://www.softwarecollections.org/>`_ for collections of software (like devtoolset 3 or python)
    * `rpmfinder <https://www.rpmfind.net/>`_ for searching rpms across mutliple systems

Setup
=====

1) Make sure that the GCC comiler is in your PATH:

::

   export PATH=/path/to/gcc:$PATH

2) Make sure that python2.7 is in your PYTHONPATH:

::

   export PYTHONPATH=/path/to/python2.7/site-packages:$PYTHONPATH

3) There are three methods of installation of ezBAMQC, from source, from setup.py, and from pypi, once prequistes are setup. 

From Source
~~~~~~~~~~

1) Download source 

2) Unpack tarball and go to the directory of the package: 

::

   tar xvfz bamqc-0.6.6.tar.gz

   cd bamqc-0.6.6

3) Run make:

::

   make

From Setup.py
~~~~~~~~~~~~

::

   python2.7 setup.py install 

From Pypi
~~~~~~~~

::

   pip2.7 install BAMqc

Usage
=====

::

   ezBAMQC [-h] -i alignment_files [alignment_files ...] -r [refgene]
   [-f [attrID]] [--rRNA [rRNA]] -o [dir] [--stranded [stranded]]
   [-q [mapq]] [-l labels [labels ...]] [-t NUMTHREADS]

optional arguments:

::

   -h, --help               show this help message and exit.
   -i, --inputFile          alignment files. Could be multiple SAM/BAM files separated by space. Required.
   -r, --refgene            gene annotation file in GTF format. Required
   -f                       the read summation at which feature level in the GTF file. DEFAULT: gene_id.
   --rRNA                   rRNA coordinates in BED format.
   -o, --outputDir          output directory. Required.
   --stranded               strandness of the library? 
                            yes : sense stranded
                            reverse : reverse stranded
                            no : not stranded
                            DEFAULT: yes.
   -q, --mapq               Minimum mapping quality (phred scaled) for an alignment to be called uniquely mapped. DEFAULT:30
   -l, --label              Labels of input files. DEFAULT:smp1 smp2 ...
   -t, --threads            Number of threads to use. DEFAULT:1

Example: 

::

   ezBAMQC -i test-data/exp_data/treat1.bam test-data/exp_data/treat2.bam test-data/exp_data/treat3.bam -r test-data/exp_data/hg9_refGene.gtf -q 30 --rRNA test-data/exp_data/hg19_rRNA.bed -o exp_output2

   Please find the example output from folder test-data.

FAQ
===
Q: Why use ezBAMQC?

A: ezBAMQC is efficient and easy to use. With one command line, it reports a comprehensive evaluation of the data with a set of plots and tables.The ability to assess multiple samples together with high efficiency make it especially useful in cases where there are a large number of samples from the same condition, genotype, or treatment. ezBAMQC was written in C++ and supports multithreading. A mouse RNA-seq sample with 120M alignments can be done in 8 minutes with 5 threads.

Q: Why the total number of reads reported by ezBAMQC does not match with samtools flagstat?

A: The difference is because of non-uniquely mapped reads or multiply aligned reads (multi-reads). Samtools flagstat counts each multiple aligment as a different reads, but ezBAMQC counts reads accoriding to the read ID, i.e., each individual read will be counted once no matter that it is a uniquely mapped read or multi-read. 

Q: What is "Low Quality Reads" ?

A: Reads marked as qc fail accoriding to SAM format or reads with mapping quality lower than the value set by the option -q will be considered as "Low Quality Reads".

Q: How the setting of option -q alter the results? 

A: Reads with low quality, i.e., did not pass -q cutoff, are only counted in Total Reads, Mapped Reads, and Mappability by mapping quality plot. The rest of the report does not include low quality reads. 

Q: Do multi-reads (non-uniquely mapped reads) have been considered in Read distribution and gene quantification?

A: No. Only uniquely mapped reads were counted. 


Acknowledgements
================

#) Samtools contributors
#) Users' valuable feedback

Copying & Distribution
======================

ezBAMQC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but *WITHOUT ANY WARRANTY*; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE*.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ezBAMQC.  If not, see `this website <http://www.gnu.org/licenses/>`_
