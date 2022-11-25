## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2015 by the deal.II authors
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

#
# Try to find the ZLIB library
#
# This module exports
#
#   ZLIB_LIBRARIES
#   ZLIB_INCLUDE_DIRS
#   ZLIB_VERSION
#

set(ZLIB_DIR "" CACHE PATH "An optional hint to a ZLIB installation")
set_if_empty(ZLIB_DIR "$ENV{ZLIB_DIR}")

if(NOT "${ZLIB_DIR}" STREQUAL "")
  set(ZLIB_ROOT ${ZLIB_DIR})
endif()
find_package(ZLIB)

set(ZLIB_VERSION ${ZLIB_VERSION_STRING})

process_feature(ZLIB
  LIBRARIES REQUIRED ZLIB_LIBRARY
  INCLUDE_DIRS REQUIRED ZLIB_INCLUDE_DIR
  CLEAR ZLIB_INCLUDE_DIR ZLIB_LIBRARY
  )
