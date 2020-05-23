Why this fork of deal.II?
=========================

This fork of deal.II is meant to add features to deal.II, the parent library of OpenFCST, in order to make the features the OpenFCST developers add to deal.II are included in the library.

The fork contains the following key branches:
    * master: This branch is to be sync with the parent project dealii/master.
    * current: This branch is sync with master, but it also contains all new deal.II features developed by the OpenFCST group. This is the version that OpenFCST uses.
    * feature branches: Each one of these branches will implement a new feature which, once it is finalized, will be merged with dealii/master via pull request. The name of the branch should be person/feature

The workflow in the fork is as follows:
a) Person A, that wants to develop a new feature for OpenFCST. Then, it would create a branch named PersonA/mesh_reader_for_GMSH from master in the OpenFCST fork (please make sure the master in the fork is up to date with dealii/master, otherwise inform the person maintaining the fork to re-synch using: git fetch origin/dealii/master).
b) Once Person A has finished the feature, submits a pull request (PR) to origin/dealii/master for review by the deal.II developers. At the same time, she/he can merge the feature into openFCST/current in our fork
c) Person A works with the deal.II developers in order to get the PR merged into deal.II.
d) Person A asks the person managing the fork to udpate OpenFCST/master with the new origin/dealii/master version with the included feature.

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

A detailed *ReadME* can be found at [./doc/readme.html](https://dealii.org/developer/readme.html)
and [./doc/users/cmake.html](https://dealii.org/developer/users/cmake.html) or at
https://www.dealii.org/.

Getting started:
----------------

The tutorial steps are located under examples/ of the installation.
Information about the tutorial steps can be found at
[./doc/doxygen/tutorial/index.html](https://dealii.org/developer/doxygen/deal.II/Tutorial.html)
or at https://www.dealii.org/.

deal.II includes support for pretty-printing deal.II objects inside GDB.
See [`contrib/utilities/dotgdbinit.py`](contrib/utilities/dotgdbinit.py) or
the new documentation page (under 'information for users') for instructions
on how to set this up.

License:
--------

Please see the file [./LICENSE.md](LICENSE.md) for details

Further information:
--------------------

For further information have a look at
[./doc/index.html](https://dealii.org/developer/index.html) or at
https://www.dealii.org.

Continuous Integration Status:
------------------------

| System | Status | More information |
| --- | --- | --- |
| Indent | [![Build Status](https://travis-ci.org/dealii/dealii.png)](https://travis-ci.org/dealii/dealii) | See https://travis-ci.org |
| Linux | [![Build Status](https://jenkins.tjhei.info/job/dealii/job/master/badge/icon)](https://jenkins.tjhei.info/job/dealii/job/master/) | See https://jenkins.tjhei.info |
| MacOS | [![Build Status](https://jenkins.tjhei.info/job/dealii-OSX/job/master/badge/icon)](https://jenkins.tjhei.info/job/dealii-OSX/job/master/) | See https://jenkins.tjhei.info |
| MacOS | [![Build Status](https://github.com/dealii/dealii/workflows/github-CI/badge.svg)](https://github.com/dealii/dealii/actions?query=workflow%3Agithub-CI) | See https://github.com/dealii/dealii/actions |
| MSVC | [![Build status](https://ci.appveyor.com/api/projects/status/e1kltrbje54ikah8/branch/master?svg=true)](https://ci.appveyor.com/project/tjhei/dealii-8th3t/branch/master) | See https://appveyor.com |
| CDash | [![cdash](https://img.shields.io/website?down_color=lightgrey&down_message=offline&label=CDash&up_color=green&up_message=up&url=https%3A%2F%2Fcdash.43-1.org%2Findex.php%3Fproject%3Ddeal.II)](https://cdash.43-1.org/index.php?project=deal.II) | Various builds and configurations on https://cdash.43-1.org/index.php?project=deal.II |

