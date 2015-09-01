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


from setuptools import setup
       

setup(name = "BAMQC",
      version = "0.5",
      description = 'Quality control tools for NGS alignment file',
      keywords='Quality control BAM file',
      packages = ['libBAMQC'],
      install_requires=['argparse','pysam>=0.8','pp'],
      scripts = ["BAMqc_mp"],
      author = "Ying Jin",
      author_email ="yjin@cshl.edu",
      license='GPLv3',
      platforms = ['Linux','MacOS'],
      url='http://hammelllab.labsites.cshl.edu/software#BAMQC',
      zip_safe = False,
      include_package_data=True )


