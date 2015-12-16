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
import sys, os, glob, fnmatch

## Added 10 Jan 2008
from distutils.core import setup, Extension
import distutils.command.install_data

## Code borrowed from wxPython's setup and config files
## Thanks to Robin Dunn for the suggestion.
## I am not 100% sure what's going on, but it works!
def opj(*args):
    path = os.path.join(*args)
    return os.path.normpath(path)

## Added 10 Jan 2008
# Specializations of some distutils command classes
class wx_smart_install_data(distutils.command.install_data.install_data):
    """need to change self.install_dir to the actual library dir"""
    def run(self):
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return distutils.command.install_data.install_data.run(self)

def find_data_files(srcdir, *wildcards, **kw):
    # get a list of all files under the srcdir matching wildcards,
    # returned in a format to be used for install_data
    def walk_helper(arg, dirname, files):
        if '.svn' in dirname:
            return
        names = []
        lst, wildcards = arg
        for wc in wildcards:
            wc_name = opj(dirname, wc)
            for f in files:
                filename = opj(dirname, f)

                if fnmatch.fnmatch(filename, wc_name) and not os.path.isdir(filename):
                    names.append(filename)
        if names:
            lst.append( (dirname, names ) )

    file_list = []
    recursive = kw.get('recursive', True)
    if recursive:
        os.path.walk(srcdir, walk_helper, (file_list, wildcards))
    else:
        walk_helper((file_list, wildcards),
                    srcdir,
                    [os.path.basename(f) for f in glob.glob(opj(srcdir, '*'))])
    return file_list

## This is a list of files to install, and where:
## Make sure the MANIFEST.in file points to all the right 
## directories too.
files = find_data_files('BAMQC/', '*.*')

from distutils.core import setup

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

###TODO HAVE TO SPLIT INTO TWO AND MAKE THE A FILE
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

BAMqc_CFLAGS = ['-fpermissive','-O3','-std=c++11','-Wno-error=declaration-after-statement'] 
BAMqc_DFLAGS = [('_FILE_OFFSET_BITS','64'),('_LARGEFILE64_SOURCE',''),('_CURSES_LIB','1')]
BAMqc_INCLUDES = ['./src/htslib']
BAMqc_HEADERS = ['./src/bamqc']
BAMqc_EXTRA = ['build/lib.linux-x86_64-2.7/htslib.so']

htslib_CFLAGS = ['-Wno-error=declaration-after-statement']
htslib_HEADERS = ['./src/htslib','./src/htslib/htslib','./src/htslib/cram']
htslib_DFLAGS = [('_FILE_OFFSET_BITS','64'),('_USE_KNETFILE','')]

setup(name = "BAMQC",
    version = "0.6.4",
    description = 'Quality control tools for NGS alignment file',
    keywords = 'Quality control BAM file',
	# make sure to add all the nessacary requires
    dependency_links=['https://gcc.gnu.org/gcc-4.8/','https://www.r-project.org/','https://cran.r-project.org/web/packages/corrplot/'],
    cmdclass = { 'install_data':    wx_smart_install_data },
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
                    extra_compile_args = htslib_CFLAGS,
                    define_macros = htslib_DFLAGS,
					libraries=['z']
                    ),
		  Extension('libBAMqc',
                    sources = BAMQC_SOURCE, 
                    extra_compile_args = BAMqc_CFLAGS,
                    include_dirs = BAMqc_HEADERS + htslib_HEADERS,
                    extra_objects = BAMqc_EXTRA,
					define_macros = BAMqc_DFLAGS
                    )
          ]        
    )
