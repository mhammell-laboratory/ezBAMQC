#!/usr/bin/env python2.7
# Setup for ezBAMQC, utilities for the Sequence Alignment/Map format.
#
#    Copyright (C) 2015 Bioinformatics Shared Resource, CSHL.
#    Portions copyright (C) 2015 Cold Spring Harbor Laboratory.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

import argparse
import sys, os, glob, fnmatch
import subprocess
from distutils.core import setup, Extension
from distutils.command.install import install as DistutilsInstall
from distutils.command.build import build as DistutilsBuild

def readme():
	with open('README.rst') as f:
		return f.read()

if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	print >> sys.stderr, "ERROR: ezBAMQC requires Python 2.7"
	sys.exit()

class Compile_Things(DistutilsBuild):
    def run(self):
        os.chdir("src/htslib/")
        subprocess.call(( "./configure"), shell=True)
        subprocess.call(( "make"), shell=True)
        os.chdir("../..")
        subprocess.call(( "make"), shell=True)
        DistutilsBuild.run(self)

class Install_Things(DistutilsInstall):
    def run(self):
        os.chdir("src/htslib/")
        subprocess.call(( "./configure"), shell=True)
        subprocess.call(( "make"), shell=True)
        os.chdir("../..")
        subprocess.call(( "make"), shell=True)
        DistutilsInstall.run(self)

BAMQC_HEADER = [
    'src/ezBAMQC/Constants.h',
    'src/ezBAMQC/Coverage_prof.h',
    'src/ezBAMQC/GeneFeatures.h',
    'src/ezBAMQC/InnerDist_prof.h',
    'src/ezBAMQC/IntervalTree.h',
    'src/ezBAMQC/Mappability.h',
    'src/ezBAMQC/parseBAM.h',
    'src/ezBAMQC/ReadDup_prof.h',
    'src/ezBAMQC/Results.h',
    'src/ezBAMQC/rRNA.h'
]

BAMQC_SOURCE = [
    'src/ezBAMQC/Coverage_prof.cpp',
    'src/ezBAMQC/GeneFeatures.cpp',
    'src/ezBAMQC/InnerDist_prof.cpp',
    'src/ezBAMQC/IntervalTree.cpp',
    'src/ezBAMQC/Mappability.cpp',
    'src/ezBAMQC/parseBAM.cpp',
    'src/ezBAMQC/ReadDup_prof.cpp',
    'src/ezBAMQC/Results.cpp',
    'src/ezBAMQC/rRNA.cpp'
]

HTSLIB_PUBLIC_HEADERS = [
	'src/htslib/bgzf.h',
	'src/htslib/faidx.h',
	'src/htslib/hfile.h',
	'src/htslib/hts.h',
	'src/htslib/hts_defs.h',
	'src/htslib/khash.h',
	'src/htslib/klist.h',
	'src/htslib/knetfile.h',
	'src/htslib/kseq.h',
	'src/htslib/ksort.h',
	'src/htslib/kstring.h',
	'src/htslib/regidx.h',
	'src/htslib/sam.h',
	'src/htslib/synced_bcf_reader.h',
	'src/htslib/tbx.h',
	'src/htslib/vcf.h',
	'src/htslib/vcf_sweep.h',
	'src/htslib/vcfutils.h'
]


HTSLIB = [
	'src/htslib/bgzf.c',
	'src/htslib/faidx.c',
	'src/htslib/hfile.c',
	'src/htslib/hfile_net.c',
	'src/htslib/hts.c',
	'src/htslib/kfunc.c',
	'src/htslib/knetfile.c',
	'src/htslib/kstring.c',
	'src/htslib/regidx.c',
	'src/htslib/sam.c',
	'src/htslib/synced_bcf_reader.c',
	'src/htslib/tbx.c',
	'src/htslib/vcf.c',
	'src/htslib/vcfutils.c',
	'src/htslib/cram/cram_codecs.c',
	'src/htslib/cram/cram_decode.c',
	'src/htslib/cram/cram_encode.c',
	'src/htslib/cram/cram_index.c',
	'src/htslib/cram/cram_io.c',
	'src/htslib/cram/cram_samtools.c',
	'src/htslib/cram/cram_stats.c',
	'src/htslib/cram/files.c',
	'src/htslib/cram/mFILE.c',
	'src/htslib/cram/md5.c',
	'src/htslib/cram/open_trace_file.c',
	'src/htslib/cram/pooled_alloc.c',
	'src/htslib/cram/rANS_static.c',
	'src/htslib/cram/sam_header.c',
	'src/htslib/cram/string_alloc.c',
	'src/htslib/cram/thread_pool.c',
	'src/htslib/cram/vlen.c',
	'src/htslib/cram/zfio.c'
]

setup(name = "ezBAMQC",
    version = "0.6.7",
    description = 'Quality control tools for NGS alignment file',
    keywords = 'Quality control BAM file',
    dependency_links=['https://gcc.gnu.org/gcc-4.8/','https://www.r-project.org/','https://cran.r-project.org/web/packages/corrplot/'],
    scripts = ["ezBAMQC"],
    author = "Ying Jin",
    author_email ="yjin@cshl.edu",
    license='GPLv3',
    platforms = ['Linux'],
    url='http://hammelllab.labsites.cshl.edu/software#BAMqc',
    long_description=readme(),
    classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: C++',
          'Operating System :: Unix',
    ],
    zip_safe = False,
    include_package_data=True,
    cmdclass={'build': Compile_Things,
              'install': Install_Things},
    package_data={'': ['libBAMqc.so'],
                  '': ['README.rst']},
    )

