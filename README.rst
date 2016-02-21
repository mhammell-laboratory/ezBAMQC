

.. image:: https://raw.githubusercontent.com/mhammell-laboratory/bamqc/master/doc/bamqc-icon.png
   :width: 200 px
   :alt: generated at codeology.braintreepayments.com/mhammell-laboratory/bamqc
   :align: right
   :target: http://codeology.braintreepayments.com/mhammell-laboratory/bamqc

BAMQC
=====

Version 0.6.5

`Github Page <https://github.com/mhammell-laboratory/bamqc>`_

`Pypi Page <https://pypi.python.org/pypi/BAMQC>`_

`MHammell Lab <http://hammelllab.labsites.cshl.edu/software>`_

Installation guide for BAMQC for from source installs
-----------------------------------------------------

`Source Code <https://github.com/mhammell-laboratory/bamqc/archive/0.6.4.tar.gz>`_

`Pypi <https://pypi.python.org/pypi?:action=display&name=BAMQC&version=0.6.5>`_

*Prerequisites:*
   * `python2.7 <https://www.python.org/download/releases/2.7/>`_
   * `R <https://www.r-project.org/>`_
   * `corrplot <https://cran.r-project.org/web/packages/corrplot/>`_
   * `GCC 4.8.1 or greater <https://gcc.gnu.org/gcc-4.8/>`_

While there are multiple methods of installing the prerequistes it may
help to look at (if using a yum based linux distro):

   * `Devtoolset-3 <https://access.redhat.com/documentation/en-US/Red_Hat_Developer_Toolset/3/html/User_Guide/sect-Red_Hat_Developer_Toolset-Install.html>`_ for GCC compilers
   * `IUS <https://ius.io/>`_ for Python2.7
   * `Software Collections <https://www.softwarecollections.org/>`_ for collections of software (like devtoolset 3 or python)
   * `rpmfinder <https://www.rpmfind.net/>`_ for searching rpms across mutliple systems

*Setup*

make sure that the GCC comiler is in your PATH

.. code:: bash
 export PATH=/path/to/gcc:$PATH


make sure that python2.7 is is in your PYTHONPATH, for permance add it to your .bashrc

.. code:: bash
 export PYTHONPATH=/path/to/python2.7/site-packages:$PYTHONPATH

There are three methods of installation of BAMqc, from source, setup.py, and from pypi, once prequistes are setup. 

*From Source*

download source 

.. code:: bash
 wget https://github.com/mhammell-laboratory/bamqc/archive/0.6.4.tar.gz
 tar xvfz bamqc-0.6.4.tar.gz
 cd bamqc-0.6.4

run make on htslib

note: make sure to have -fPIC in your cflags ie: ``export CFLAGS="$CFLAGS -fPIC"``

.. code:: bash
 
 cd src/htslib
 make 

run make 

.. code:: bash

 cd ../..
 make

*From Setup.py*

.. code:: bash

 python2.7 setup.py install 

*From Pypi*

.. code:: bash

 pip2.7 install BAMqc

Contacts
--------

Ying Jin: yjin@cshl.edu

Acknowledgements goes to
------------------------

#) Samtools and pysam contributors
#) Molly Hammell and members of her laboratory at Cold Spring Harbor Laboratory
#) Users' valuable feedback

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

