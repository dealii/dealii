## ---------------------------------------------------------------------
##
## Copyright (C) 2014 by the deal.II authors
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

#
# Houston, we have a problem: CMake ships its own FindBZip2.cmake module.
# Unfortunately we want to call DEAL_II_PACKAGE_HANDLE. Therefore, use the
# original find module and do a dummy call to DEAL_II_PACKAGE_HANDLE:
#

LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
FIND_PACKAGE(BZip2)
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
