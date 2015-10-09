import sys

"""
Setup script for BAMQC  -- Comprehensive QC package for NGS data alignment file
Copyright (c) 2015 Ying Jin <yjin@cshl.edu>
This code is free software;you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).
"""

if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	print >> sys.stderr, "ERROR: RSeQC requires Python 2.7"
	sys.exit()

from setuptools import setup, Extension
from Cython.Build import cythonize

command_classes = {}

try:
    import Cython.Distutils
    command_classes['build_ext'] = Cython.Distutils.build_ext
except:
    pass

setup(name = "BAMQC",
      version = "0.5.1",
      description = 'Quality control tools for NGS alignment file',
      #long_description=readme(),
      keywords='Quality control BAM file',
      packages = ['libBAMQC'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Natural Language :: English',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
		  'Intended Audience :: Science/Research',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: C++',
          'Operating System :: MacOS',
          'Operating System :: Unix'
      ],
      install_requires=['argparse','pysam>=0.8'],
      scripts = ["BAMqc_mp"],
      author = "Ying Jin",
      author_email ="yjin@cshl.edu",
      license='GPLv3',
      platforms = ['Linux','MacOS'],
      url='http://hammelllab.labsites.cshl.edu/software#BAMQC',
      zip_safe = False,
      include_package_data=True,
      cmdclass=command_classes,
      ext_modules = cythonize([Extension( "PyGeneFeatures",
                              sources = [ "src/PyGeneFeatures.cpp","src/IntervalTree.cpp","src/GeneFeatures.cpp", "src/rRNA.cpp"],
                              language="c++",
                              extra_compile_args=["-std=c++11"],
                              extra_link_args=["-std=c++11"],
                              include_dirs=["src"])]),
      )
