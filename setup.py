#!/usr/bin/env python2.7
# Setup for BAMqc, utilities for the Sequence Alignment/Map format.
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
import sys, os, subprocess
from distutils.core import setup, Extension

def readme():
	with open('README.rst') as f:
		return f.read()

if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	print >> sys.stderr, "ERROR: BAMQC requires Python 2.7"
	sys.exit()

BAMQC_HEADER = [
    'src/bamqc/Constants.h',
    'src/bamqc/Coverage_prof.h',
    'src/bamqc/GeneFeatures.h',
    'src/bamqc/InnerDist_prof.h',
    'src/bamqc/IntervalTree.h',
    'src/bamqc/Mappability.h',
    'src/bamqc/parseBAM.h',
    'src/bamqc/ReadDup_prof.h',
    'src/bamqc/Results.h',
    'src/bamqc/rRNA.h'
]

BAMQC_SOURCE = [
    'src/bamqc/Coverage_prof.cpp',
    'src/bamqc/GeneFeatures.cpp',
    'src/bamqc/InnerDist_prof.cpp',
    'src/bamqc/IntervalTree.cpp',
    'src/bamqc/Mappability.cpp',
    'src/bamqc/parseBAM.cpp',
    'src/bamqc/ReadDup_prof.cpp',
    'src/bamqc/Results.cpp',
    'src/bamqc/rRNA.cpp'
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

BAMqc_CFLAGS = ['-fpermissive','-O3','-std=c++11'] 
BAMqc_DFLAGS = ['-D_FILE_OFFSET_BITS=64','-D_LARGEFILE64_SOURCE','-D_CURSES_LIB=1']
BAMqc_INCLUDES = ['./src/htslib']
BAMqc_HEADERS = ['./src/bamqc']

htslib_HEADERS = ['./src/htslib','./src/htslib/htslib','./src/htslib/cram']

setup(name = "BAMQC",
    version = "0.6.0",
    description = 'Quality control tools for NGS alignment file',
    keywords = 'Quality control BAM file',
	# make sure to add all the nessacary requires
    dependency_links=['https://gcc.gnu.org/gcc-4.8/','https://www.r-project.org/','https://cran.r-project.org/web/packages/corrplot/'],
    scripts = ["BAMqc"],
    author = "Ying Jin",
    author_email ="yjin@cshl.edu",
    license='GPLv3',
    platforms = ['Linux'],
    url='http://hammelllab.labsites.cshl.edu/software#BAMqc',
    long_description=readme(),
    classifiers=[
          'Development Status :: 4 - Beta',
          'Natural Language :: English',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: C++',
          'Operating System :: Unix',
    ],
    zip_safe = False,
    include_package_data=True,
    ext_modules = [ 
	      Extension('htslib',
                    sources = HTSLIB,
                    include_dirs = htslib_HEADERS,
                    language = 'c++'
                    ),
		  Extension('libBAMqc',
                    sources = BAMQC_SOURCE, 
                    extra_compile_args = BAMqc_CFLAGS + BAMqc_DFLAGS,
                    include_dirs = BAMqc_HEADERS + htslib_HEADERS,
                    #this needs to be made relative or take into acount final destination, testing is needed:
                    extra_objects = ['build/lib.linux-x86_64-2.7/htslib.so'],
                    language = 'c++',
                    )
          ]        
    )
