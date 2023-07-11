[![workflows/indent](https://github.com/dealii/dealii/workflows/indent/badge.svg)](https://github.com/dealii/dealii/actions?query=workflow%3Aindent)
[![workflows/tidy](https://github.com/dealii/dealii/workflows/tidy/badge.svg)](https://github.com/dealii/dealii/actions?query=workflow%3Atidy)
[![workflows/github-linux](https://github.com/dealii/dealii/workflows/github-linux/badge.svg)](https://github.com/dealii/dealii/actions?query=workflow%3Agithub-linux)
[![workflows/github-OSX](https://github.com/dealii/dealii/workflows/github-OSX/badge.svg)](https://github.com/dealii/dealii/actions?query=workflow%3Agithub-OSX)
[![workflows/github-windows](https://github.com/dealii/dealii/workflows/github-windows/badge.svg)](https://github.com/dealii/dealii/actions?query=workflow%3Agithub-windows)
[![workflows/github-docker](https://github.com/dealii/dealii/workflows/github-docker/badge.svg)](https://github.com/dealii/dealii/actions?query=workflow%3Agithub-docker)
[![jenkins/dealii-serial](https://ci.tjhei.info/job/dealii-serial/job/master/badge/icon?style=plastic&subject=jenkins-serial)](https://ci.tjhei.info/job/dealii-serial/job/master/)
[![jenkins/dealii-mpi](https://ci.tjhei.info/job/dealii-mpi/job/master/badge/icon?style=plastic&subject=jenkins-MPI)](https://ci.tjhei.info/job/dealii-mpi/job/master/)
[![jenkins/dealii-osx](https://ci.tjhei.info/job/dealii-osx/job/master/badge/icon?style=plastic&subject=jenkins-OSX)](https://ci.tjhei.info/job/dealii-osx/job/master/)
[![jenkins/dealii-ampere](https://ci.tjhei.info/job/dealii-ampere/job/master/badge/icon?style=plastic&subject=jenkins-ampere)](https://ci.tjhei.info/job/dealii-ampere/job/master/)

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

A detailed *ReadME* can be found at [./doc/readme.html](https://dealii.org/developer/readme.html),
[./doc/users/cmake_user.html](https://dealii.org/developer/users/cmake_user.html),
or at https://www.dealii.org/.

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

Docker Images:
-------------

Docker images based on the Ubuntu operating system are available on
[Docker Hub](https://hub.docker.com/repository/docker/dealii/dealii). You can 
use any of the available version 
([list of available tags](https://hub.docker.com/repository/docker/dealii/dealii/tags)) 
by running, for example:

    $ docker run --rm -t -i dealii/dealii:master-focal

The above command would drop you into an isolated environment, in which you 
will find the latest version of deal.II (master development branch) installed
under `/usr/local`.

