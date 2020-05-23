## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

########################################################################
#                                                                      #
#                   Check for various compiler bugs:                   #
#                                                                      #
########################################################################


#
# In intel (at least 13.1 and 14), vectorization causes
# wrong code. See https://code.google.com/p/dealii/issues/detail?id=156
# or tests/hp/solution_transfer.cc
# A work-around is to disable all vectorization.
#
# - Timo Heister, 2013, 2015
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0.3" )
  ENABLE_IF_SUPPORTED(DEAL_II_CXX_FLAGS_RELEASE "-no-vec")
ENDIF()


#
# Intel 16.0.1 produces wrong code that creates a race condition in
# tests/fe/curl_curl_01.debug but 16.0.2 is known to work. Blacklist this
# version. Also see github.com/dealii/dealii/issues/2203
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "Intel" AND CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "16.0.1" )
  MESSAGE(FATAL_ERROR "Intel compiler version 16.0.1 is not supported, please update to 16.0.2 or newer!")
ENDIF()
