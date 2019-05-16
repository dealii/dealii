## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
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
# Try to find Mesquite
#
# This module exports
#
#   MESQUITE_INCLUDE_DIR
#   MESQUITE_LIBRARY
#   MESQUITE_VERSION
#   MESQUITE_VERSION_MAJOR
#   MESQUITE_VERSION_MINOR
#   MESQUITE_VERSION_SUBMINOR
#

SET(MESQUITE_DIR "" CACHE PATH "An optional hint to a Mesquite installation")
SET_IF_EMPTY(MESQUITE_DIR "$ENV{MESQUITE_DIR}")

DEAL_II_FIND_PATH(MESQUITE_INCLUDE_DIR
  NAMES Mesquite_all_headers.hpp
  HINTS ${MESQUITE_DIR}
  PATH_SUFFIXES include
  )

DEAL_II_FIND_LIBRARY(MESQUITE_LIBRARY
  NAMES mesquite
  HINTS ${MESQUITE_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_LIBRARY(MESQUITE_UTILITY_LIBRARY
  NAMES msqutil
  HINTS ${MESQUITE_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# Look for mesquite_version.h - we'll query it to determine the library version:
#

DEAL_II_FIND_FILE(MESQUITE_VERSION_H mesquite_version.h
  HINTS ${MESQUITE_INCLUDE_DIR}
  NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
  )


IF(EXISTS ${MESQUITE_VERSION_H})
  FILE(STRINGS "${MESQUITE_VERSION_H}" MESQUITE_VERSION_MAJOR_STRING
    REGEX "#define MSQ_VERSION_MAJOR"
    )
  FILE(STRINGS "${MESQUITE_VERSION_H}" MESQUITE_VERSION_MINOR_STRING
    REGEX "#define MSQ_VERSION_MINOR"
    )
  FILE(STRINGS "${MESQUITE_VERSION_H}" MESQUITE_VERSION_PATCH_STRING
    REGEX "#define MESQUITE_VERSION_PATCH"
    )

  STRING(REGEX REPLACE "^.*MSQ_VERSION_MAJOR[ ]+([0-9]+).*" "\\1"
    MESQUITE_VERSION_MAJOR "${MESQUITE_VERSION_MAJOR_STRING}"
    )
  STRING(REGEX REPLACE "^.*MSQ_VERSION_MINOR[ ]+([0-9]+).*" "\\1"
    MESQUITE_VERSION_MINOR "${MESQUITE_VERSION_MINOR_STRING}"
    )
  STRING(REGEX REPLACE "^.*MSQ_VERSION_MAJOR[ ]+([0-9]+).*" "\\1"
    MESQUITE_VERSION_PATCH "${MESQUITE_VERSION_PATCH_STRING}"
    )

  IF(NOT ${MESQUITE_VERSION_PATCH} STREQUAL "")
    SET(MESQUITE_VERSION
      "${MESQUITE_VERSION_MAJOR}.${MESQUITE_VERSION_MINOR}.${MESQUITE_VERSION_PATCH}"
      )
  ELSE()
    SET(MESQUITE_VERSION
      "${MESQUITE_VERSION_MAJOR}.${MESQUITE_VERSION_MINOR}"
      )
  ENDIF()
ENDIF()


DEAL_II_PACKAGE_HANDLE(MESQUITE
  LIBRARIES 
    REQUIRED MESQUITE_LIBRARY
    OPTIONAL MESQUITE_UTILITY_LIBRARY
  INCLUDE_DIRS REQUIRED MESQUITE_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED MESQUITE_INCLUDE_DIR
  CLEAR MESQUITE_INCLUDE_DIR MESQUITE_LIBRARY MESQUITE_UTILITY_LIBRARY
  )
