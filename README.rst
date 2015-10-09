BAMQC
=====

`Github Page <https://github.com/mhammell-laboratory/bamqc>`_

`Pypi Page <https://pypi.python.org/pypi/BAMQC>`_

`MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

Installation guide for BAMQC for from source installs
-----------------------------------------------------

*Prerequisites:*
   * python2.7
   * GCC 4.9 or greater

Below is an example of installing BAMQC on Linux system using BASH. You need to change '--root' directory, PYTHONPATH and PATH accordingly

::
     tar zxf BAMQC-VERSION.tar.gz
     cd BAMQC-VERSION
     python setup.py install will install BAMQC in system level. require root previledge
     python setup.py install --prefix=/home/user/BAMQC     #will install BAMQC at user specified location
     export PYTHONPATH=/home/user/BAMQC/usr/local/lib/python2.7/site-packages:$PYTHONPATH.     #setup PYTHONPATH, so that BAMQC knows where to import modules
     export PATH=/home/user/BAMQC/usr/local/bin:$PATH     #setup PATH, so that system knows where to find executable file 


*NOTE:*

* To install BAMQC on MAC OSX, user need to download and install Xcode beforehand.
* To produce graphical outputs, R must be installed.
* If the installation failed with error like: /usr/bin/ld: cannot find -lz, you may need to install a shared zlib library on your system. 

Contacts
--------

Ying Jin: yjin@cshl.edu

Acknowledgements goes to
------------------------

1. Samtools and pysam contributors
2. Molly Hammell and members in Molly Hammell's Laboratory (Cold Spring Harbor Laboratory)
3. Users' valuable feedback

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
along with TEToolKit.  If not, see `this website <http://www.gnu.org/licenses/>`_.