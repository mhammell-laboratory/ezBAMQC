BAMqc
======

Version: X.X.X

Quality analysis of sequencing data using aligned files (BAM)

`MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

Created by Ying Jin & Molly Hammell, August 2014

Copyright (C) 2014-2015 Ying Jin & Molly Hammell

Contact: Ying Jin (yjin@cshl.edu)


Requirements
------------

Python:      2.6.x or 2.7.x (not tested in Python 3.x)


Installation
-------------

1. Download compressed tarball.
2. Unpack tarball.
3. Navigate into unpacked directory.
4. Run the following::

    $ python setup.py install

If you want to install locally (e.g. /local/home/usr),
run this command instead::

    $ python setup.py install --prefix /local/home/usr

*NOTE* In the above example, you must add

    /local/home/usr/bin

to the PATH variable, and

     /local/home/usr/lib/python2.X/site-packages 

to the PYTHONPATH variable, where python2.X refers to the 
python version (e.g. python2.7 if using python version 2.7.x).


Usage
-----

    usage: BAMqc [-h] -f alignment_files [alignment_files ...] -r [refgene] -o
             [dir] [-i [transript_Index]] [-q [mapq]] [-l [lb]] [-u [ub]]
             [-s [stepsize]] [-t labels [labels ...]]

    Optional arguments:
      -h, --help            show this help message and exit
      -f alignment_files [alignment_files ...], --inputFile alignment_files [alignment_files ...]
                        Alignment files. Could be multiple BAM files separated
                        by space.
      -r [refgene], --refgene [refgene]
                        refGene BED12 file.
      -o [dir], --outputDir [dir]
                        output directory.
      -i [transript_Index], --index [transript_Index]
                        Transcriptome index file.
      -q [mapq], --mapq [mapq]
                        Minimum mapping quality (phred scaled) for an
                        alignment to be called uniquely mapped. DEFAULT:30
      -l [lb], --lowBound [lb]
                        Lower bound for plotting insert size distribution.
                        DEFAULT:-250
      -u [ub], --upperBound [ub]
                        Upper bound for plotting insert size distribution.
                        DEFAULT:250
      -s [stepsize], --stepSize [stepsize]
                        Upper bound for plotting insert size distribution.
                        DEFAULT:5
      -t labels [labels ...], --label labels [labels ...]
                        Labels of input files. DEFAULT:smp1 smp2 ...


Example Command Lines
---------------------

    BAMqc -f treat1.bam treat2.bam treat3.bam -r mm9_refGene.bed

