import sys
from setuptools import setup, Extension
import distutils.ccompiler

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

HTSLIB = [
    'htslib/',
	'htslib/'
]

command_classes = {}

setup(name = "BAMQC",
    version = "0.6.0",
    description = 'Quality control tools for NGS alignment file',
    keywords='Quality control BAM file',
    packages = ['BAMqc'],
    install_requires=['argparse','pysam>=0.8'],
    scripts = ["BAMqc"],
    author = "Ying Jin",
    author_email ="yjin@cshl.edu",
    license='GPLv3',
    platforms = ['Linux','MacOS'],
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
