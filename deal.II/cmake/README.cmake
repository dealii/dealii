## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------
This folder contains the deal.II CMake build system

Extensive documentation can be found at /doc/development/cmake-internals.html

It is structured as follows:

./checks
========

Contains checks for platform features and compiler bugs and features

./config
=======

Contains configuration templates for
  - the project configuration (deal.IIConfig.cmake)
  - the legacy Make.global_options mechanism
  - the C++ template expansion mechanism (template-arguments)

./configure
===========

Contains files configure_<feature>.cmake for configuration and setup of
all features the deal.II library supports

./macros
========

CMake script macros for several purposes

./modules
=========

Contains Find<Library>.cmake modules for finding external libraries

./scripts
=========

Contains script files needed for the build system, notably expand_instantiations


./setup_*.cmake
===============

Setup files included by the top level CMakeLists.txt file
