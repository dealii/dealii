## ---------------------------------------------------------------------
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


#
# CMakes Ninja generator is currently incompatible with everything but gcc
# and clang.
#
# - Matthias Maier, 2013
#

IF( CMAKE_GENERATOR MATCHES "Ninja" AND NOT
    ( CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
      CMAKE_CXX_COMPILER_ID MATCHES "Clang"  ) )
  MESSAGE(FATAL_ERROR "\n"
    "Error!\n"
    "The CMAKE_GENERATOR \"${CMAKE_GENERATOR}\" "
    "currently only supports the GNU and Clang C++ compilers, but "
    "\"${CMAKE_CXX_COMPILER_ID}\" was found.\n\n"
    )
ENDIF()
