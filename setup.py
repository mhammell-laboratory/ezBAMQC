import sys
from setuptools import setup, Extension
import distutils.ccompiler

###THIS FEATURE DOESN'T work yet###

"""
Setup script for BAMQC  -- Comprehensive QC package for NGS data alignment file
Copyright (c) 2015 Ying Jin <yjin@cshl.edu>
This code is free software;you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).
"""
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

BAMQC_DOCS = [
    'doc/CONTACTS',
    'doc/COPYING',
    'doc/INSTALL',
    'doc/THANKS'
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
	'src/bgzf.c',
	'src/faidx.c',
	'src/hfile_internal.h',
	'src/hfile.c',
	'src/hfile_irods.c',
	'src/hfile_net.c',
	'src/hts.c',
	'src/knetfile.c',
	'src/kstring.c',
	'src/regidx.c',
	'src/sam.c',
	'src/synced_bcf_reader.c',
	'src/tbx.c',
	'src/vcf.c',
	'src/vcf_sweep.c',
	'src/vcfutils.c',
	'src/cram/cram.h',
	'src/cram/cram_codecs.c',
	'src/cram/cram_codecs.h',
	'src/cram/cram_decode.c',
	'src/cram/cram_decode.h',
	'src/cram/cram_encode.c',
	'src/cram/cram_encode.h',
	'src/cram/cram_index.c',
	'src/cram/cram_index.h',
	'src/cram/cram_io.c',
	'src/cram/cram_io.h',
	'src/cram/cram_samtools.c',
	'src/cram/cram_samtools.h',
	'src/cram/cram_stats.c',
	'src/cram/cram_stats.h',
	'src/cram/cram_structs.h',
	'src/cram/files.c',
	'src/cram/mFILE.c',
	'src/cram/mFILE.h',
	'src/cram/md5.c',
	'src/cram/md5.h',
	'src/cram/misc.h',
	'src/cram/open_trace_file.c',
	'src/cram/open_trace_file.h',
	'src/cram/os.h',
	'src/cram/pooled_alloc.c',
	'src/cram/pooled_alloc.h',
	'src/cram/sam_header.c',
	'src/cram/sam_header.h',
	'src/cram/string_alloc.c',
	'src/cram/string_alloc.h',
	'src/cram/thread_pool.c',
	'src/cram/thread_pool.h',
	'src/cram/vlen.c',
	'src/cram/vlen.h',
	'src/cram/zfio.c',
	'src/cram/zfio.h'
]

command_classes = {}

setup(name = "BAMQC",
    version = "0.6.0",
    description = 'Quality control tools for NGS alignment file',
    keywords='Quality control BAM file',
    packages = ['BAMqc'],
    install_requires=['argparse',''],
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
    cmdclass=command_classes,
    )
