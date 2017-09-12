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
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Try to find the BZIP2 library
#
# This module exports
#
#   BZIP2_LIBRARIES
#   BZIP2_INCLUDE_DIRS
#   BZIP2_VERSION
#

SET(BZIP2_DIR "" CACHE PATH "An optional hint to a BZIP2 installation")
SET_IF_EMPTY(BZIP2_DIR "$ENV{BZIP2_DIR}")

SET(_cmake_prefix_path_backup "${CMAKE_PREFIX_PATH}")

# temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

SET(CMAKE_PREFIX_PATH ${BZIP2_DIR} ${_cmake_prefix_path_backup})

FIND_PACKAGE(BZip2)

SET(CMAKE_PREFIX_PATH ${_cmake_prefix_path_backup})
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

SET(BZIP2_VERSION ${BZIP2_VERSION_STRING})
SET(_bzip2_libraries ${BZIP2_LIBRARIES})

DEAL_II_PACKAGE_HANDLE(BZIP2
  LIBRARIES REQUIRED _bzip2_libraries
  INCLUDE_DIRS REQUIRED BZIP2_INCLUDE_DIR
  CLEAR
    BZIP2_INCLUDE_DIR BZIP2_LIBRARY_DEBUG BZIP2_LIBRARY_RELEASE
    BZIP2_NEED_PREFIX
  )
