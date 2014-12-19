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
# Try to find the ZLIB library
#
# This module exports
#
#   ZLIB_LIBRARIES
#   ZLIB_INCLUDE_DIRS
#   ZLIB_VERSION
#

SET(ZLIB_DIR "" CACHE PATH "An optional hint to a ZLIB installation")
SET_IF_EMPTY(ZLIB_DIR "$ENV{ZLIB_DIR}")

#
# Houston, we have a problem: CMake ships its own FindZLIB.cmake module.
# Unfortunately we want to call DEAL_II_PACKAGE_HANDLE. Therefore, use the
# original find module and do a dummy call to DEAL_II_PACKAGE_HANDLE:
#

IF(NOT "${ZLIB_DIR}" STREQUAL "")
  SET(ZLIB_ROOT ${ZLIB_DIR})
ENDIF()
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
FIND_PACKAGE(ZLIB)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

SET(ZLIB_VERSION ${ZLIB_VERSION_STRING})

DEAL_II_PACKAGE_HANDLE(ZLIB
  LIBRARIES REQUIRED ZLIB_LIBRARY
  INCLUDE_DIRS REQUIRED ZLIB_INCLUDE_DIR
  CLEAR ZLIB_INCLUDE_DIR ZLIB_LIBRARY
  )
