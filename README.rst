

.. image:: https://raw.githubusercontent.com/mhammell-laboratory/bamqc/master/doc/bamqc-icon.png
   :width: 200 px
   :alt: generated at codeology.braintreepayments.com/mhammell-laboratory/bamqc
   :align: right
   :target: http://codeology.braintreepayments.com/mhammell-laboratory/bamqc



BAMQC
=====

Version 0.6.4

BAMQC is a tool to check the quality of either one or many mapped next-generation-sequencing
datasets. It conducts comprehensive evaluations of aligned sequencing data from multiple aspects including: clipping
profile, mapping quality distribution, mapped read length distribution, genomic/transcriptomic mapping distribution, inner
distance distribution (for paired-end reads), ribosomal RNA contamination, transcript 5’ and 3’ end bias, transcription
dropout rate, sample correlations, sample reproducibility, sample variations. It outputs a set of tables and plots and one HTML
page that contains a summary of the results. Many metrics are designed for RNA-seq data specifically, but BAMQC can be
applied to any mapped sequencing dataset such as RNA-seq, CLIP-seq, GRO-seq, ChIP-seq, DNA-seq and so on.


`Github Page <https://github.com/mhammell-laboratory/bamqc>`_

`Pypi Page <https://pypi.python.org/pypi/BAMQC>`_

`MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

Created by Ying Jin, David Molik, and Molly Hammell, 2015

Contact: Ying Jin (yjin@cshl.edu)

Installation guide for BAMQC for from source installs
-----------------------------------------------------

`Source Code <https://github.com/mhammell-laboratory/bamqc/archive/0.6.4.tar.gz>`_

`Pypi <https://pypi.python.org/pypi?:action=display&name=BAMQC&version=0.6.4>`_

*Prerequisites:*
   * `python2.7 <https://www.python.org/download/releases/2.7/>`_
   * `R <https://www.r-project.org/>`_
   * `corrplot <https://cran.r-project.org/web/packages/corrplot/>`_
   * `GCC 4.8.1 or greater <https://gcc.gnu.org/gcc-4.8/>`_

*While there are multiple methods of installing the prerequistes it may help to look at (if using a yum based linux distro):*
   * `Devtoolset-3 <https://access.redhat.com/documentation/en-US/Red_Hat_Developer_Toolset/3/html/User_Guide/sect-Red_Hat_Developer_Toolset-Install.html>`_ for GCC compilers
   * `IUS <https://ius.io/>`_ for Python2.7
   * `Software Collections <https://www.softwarecollections.org/>`_ for collections of software (like devtoolset 3 or python)
   * `rpmfinder <https://www.rpmfind.net/>`_ for searching rpms across mutliple systems

*Setup*

1) Make sure that the GCC comiler is in your PATH:
  .. code:: bash

   export PATH=/path/to/gcc:$PATH

2) Make sure that python2.7 is in your PYTHONPATH:
  .. code:: bash

   export PYTHONPATH=/path/to/python2.7/site-packages:$PYTHONPATH

3) There are three methods of installation of BAMQC, from source, setup.py, and from pypi, once prequistes are setup. 
 -*From Source*
  1) Download source 
  2) Unpack tarball and go to the directory of the package: 
   .. code:: bash

    tar xvfz bamqc-0.6.4.tar.gz

    cd bamqc-0.6.4
  3) Run make:
   .. code:: bash

    make
 -*From Setup.py*
  .. code:: bash

   python2.7 setup.py install 
 -*From Pypi*
  .. code:: bash

   pip2.7 install BAMqc

*Usage*

 **BAMQC [-h] -i alignment_files [alignment_files ...] -r [refgene]**

              **[-f [attrID]] [--rRNA [rRNA]] -o [dir] [--stranded [stranded]]**

              **[-q [mapq]] [-l labels [labels ...]] [-t NUMTHREADS]**

 optional arguments:

  **-h, --help**            show this help message and exit

  -i alignment_files [alignment_files ...], --inputFile alignment_files [alignment_files ...]
                        Alignment files. Could be multiple SAM/BAM files
                        separated by space. Required.

  -r [refgene], --refgene [refgene] gene annotation file in GTF format. Required

  -f [attrID]           The read summation at which feature level in the GTF
                        file. DEFAULT: gene_id.

  --rRNA [rRNA]         rRNA coordinates in BED format.

  -o [dir], --outputDir [dir] output directory. Required.

  --stranded [stranded] strandness of the library? 

                        yes : sense stranded

                        reverse : reverse stranded

                        no : not stranded

                        DEFAULT: yes.

  -q [mapq], --mapq [mapq] Minimum mapping quality (phred scaled) for an alignment to be called uniquely mapped. DEFAULT:30

  -l labels [labels ...], --label labels [labels ...] Labels of input files. DEFAULT:smp1 smp2 ...

  -t NUMTHREADS, --threads NUMTHREADS Number of threads to use. DEFAULT:1

 Example: BAMQC -i treat1.bam treat2.bam treat3.bam -r mm9_refGene.gtf -q 30 --rRNA mm9_rRNA.bed -o bamqc_out


Acknowledgements goes to
------------------------

#) Samtools and pysam contributors
#) Users' valuable feedback

Copying & distribution
----------------------

BAMQC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but *WITHOUT ANY WARRANTY*; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE*.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BAMQC.  If not, see `this website <http://www.gnu.org/licenses/>`_.

