What is deal.II?
================

deal.II is a C++ program library targeted at the computational solution
of partial differential equations using adaptive finite elements. It uses
state-of-the-art programming techniques to offer you a modern interface
to the complex data structures and algorithms required.

For the impatient:
------------------

Let's say you've unpacked the .tar.gz file into a directory /path/to/dealii/sources. 
Then configure, compile, and install the deal.II library with:

    $ mkdir build
    $ cd build
    $ cmake -DCMAKE_INSTALL_PREFIX=/path/where/dealii/should/be/installed/to /path/to/dealii/sources
    $ make install    (alternatively $ make -j<N> install)
    $ make test

To build from the repository, execute the following commands first:

    $ git clone https://github.com/dealii/dealii
    $ cd dealii

Then continue as before.

A detailed *ReadME* can be found at ./doc/readme.html and
./doc/users/cmake.html or at http://www.dealii.org/.

Getting started:
----------------

The tutorial steps are located under examples/ of the installation.
Information about the tutorial steps can be found at
doc/doxygen/tutorial/index.html or at http://www.dealii.org/.

deal.II includes support for pretty-printing deal.II objects inside GDB. See
`contrib/utilities/dotgdbinit.py` or the new documentation page (under
'information for users') for instructions on how to set this up.

License:
--------

Please see the file ./LICENSE for details

Further information:
--------------------

For further information have a look at ./doc/index.html or at
http://www.dealii.org.

Continuous Integration Status:
------------------------

[![Build Status](https://travis-ci.org/dealii/dealii.png)](https://travis-ci.org/dealii/dealii)
