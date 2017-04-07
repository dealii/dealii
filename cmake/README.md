This folder contains the deal.II CMake build system
===================================================

Extensive documentation can be found at
`/doc/developers/cmake-internals.html`

It is structured as follows:

checks
------

Contains checks for platform features and compiler bugs and features

config
------

Contains configuration templates for

  * the project configuration (deal.IIConfig.cmake)
  * the C++ template expansion mechanism (template-arguments)

configure
---------

Contains files `configure_<feature>.cmake` for configuration and setup of
all features the deal.II library supports

macros
------

CMake script macros for several purposes

modules
-------

Contains `Find<Library>.cmake` modules for finding external libraries

scripts
-------

Contains script files needed for the build system, notably expand_instantiations

setup_*.cmake
-------------

Setup files included by the top level `CMakeLists.txt` file
