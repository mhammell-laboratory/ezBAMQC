BAMQC
=====

*Version 0.5.4*

`Github Page <https://github.com/mhammell-laboratory/bamqc>`_

`Pypi Page <https://pypi.python.org/pypi/BAMQC>`_

`MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

Installation guide for BAMQC for from source installs
-----------------------------------------------------

`Source Code <https://github.com/mhammell-laboratory/bamqc/archive/0.5.4.tar.gz>`_

`With Prebuilt Binary <https://github.com/mhammell-laboratory/bamqc/releases/download/0.5.4/BAMQC-0.5.4.tar.gz>`_

`Pypi <https://pypi.python.org/pypi/BAMQC>`_

*Prerequisites:*
   * python2.7
   * R (corrplot package)
   * GCC 4.8.1 or greater (Linux), Xcode 4.2 or greater (MacOSX)

Below is an example of installing BAMQC on Linux system using BASH. You need to change '--prefix' directory, PYTHONPATH and PATH accordingly

::

    tar zxf BAMQC-VERSION.tar.gz
    cd BAMQC-VERSION

To install BAMQC at the system level (which will require root privileges):

::

    python setup.py install

To install BAMQC in a custom location (e.g. /home/user/BAMQC):

::

    export PYTHONPATH=/home/user/BAMQC/lib/python2.7/site-packages:$PYTHONPATH.
    # This sets up PYTHONPATH so that BAMQC knows where to install, later import, BAMQC modules.
    export PATH=/home/user/BAMQC/bin:$PATH
    # This sets up PATH, so that system knows where to find the executable file.
    python setup.py install --prefix=/home/user/BAMQC
    # This will install BAMQC at the user specified location (/home/user/BAMQC)


*NOTE:*

* To produce graphical outputs, R (and the corrplot R package) must be installed.
* If the installation failed with error like: /usr/bin/ld: cannot find -lz, you may need to install a shared zlib library on your system.

Contacts
--------

Ying Jin: yjin@cshl.edu

Acknowledgements goes to
------------------------

1. Samtools and pysam contributors
2. Molly Hammell and members of her laboratory at Cold Spring Harbor Laboratory
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
along with BAMQC.  If not, see `this website <http://www.gnu.org/licenses/>`_.

