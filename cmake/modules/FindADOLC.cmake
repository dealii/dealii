## ---------------------------------------------------------------------
##
## Copyright (C) 2016 - 2018 by the deal.II authors
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

# Note: This file is adapted from https://raw.githubusercontent.com/joaoleal/CppADCodeGen/master/cmake/FindADOLC.cmake

#
# Try to find Adolc
#
# This module exports
#
#   ADOLC_INCLUDE_DIR
#   ADOLC_LIBRARY
#   ADOLC_WITH_ADVANCED_BRANCHING
#   ADOLC_WITH_ATRIG_ERF
#   ADOLC_WITH_BOOST_ALLOCATOR
#

SET(ADOLC_DIR "" CACHE PATH "An optional hint to an ADOL-C installation")
SET_IF_EMPTY(ADOLC_DIR "$ENV{ADOLC_DIR}")

DEAL_II_FIND_PATH(ADOLC_INCLUDE_DIR
  NAMES adolc/adolc.h
  HINTS ${ADOLC_DIR}
  PATH_SUFFIXES include
  )

DEAL_II_FIND_LIBRARY(ADOLC_LIBRARY
  NAMES adolc
  HINTS ${ADOLC_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# Look for adolc_settings.h - we'll query it to determine supported features:
#

DEAL_II_FIND_FILE(ADOLC_SETTINGS_H adolc_settings.h
  HINTS ${ADOLC_INCLUDE_DIR} "${ADOLC_INCLUDE_DIR}/adolc/internal"
  NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
  )

IF(EXISTS ${ADOLC_SETTINGS_H})
  #
  # Check whether ADOL-C is configured with extra trig functions
  #
  FILE(STRINGS "${ADOLC_SETTINGS_H}" ADOLC_ATRIG_ERF_STRING
    REGEX "#define ATRIG_ERF"
    )
  IF(NOT "${ADOLC_ATRIG_ERF_STRING}" STREQUAL "")
    SET(ADOLC_WITH_ATRIG_ERF TRUE)
  ELSE()
    SET(ADOLC_WITH_ATRIG_ERF FALSE)
  ENDIF()

  #
  # Check whether ADOL-C is configured with advanced branching
  #
  FILE(STRINGS "${ADOLC_SETTINGS_H}" ADOLC_ADVANCED_BRANCHING_STRING
    REGEX "#define ADOLC_ADVANCED_BRANCHING"
    )
  IF(NOT "${ADOLC_ADVANCED_BRANCHING_STRING}" STREQUAL "")
    SET(ADOLC_WITH_ADVANCED_BRANCHING TRUE)
  ELSE()
    SET(ADOLC_WITH_ADVANCED_BRANCHING FALSE)
  ENDIF()

  #
  # Check whether ADOL-C is configured to use the Boost pool allocator
  #
  FILE(STRINGS "${ADOLC_SETTINGS_H}" ADOLC_BOOST_POOL_STRING
    REGEX "#define USE_BOOST_POOL 1"
    )
  IF(NOT "${ADOLC_BOOST_POOL_STRING}" STREQUAL "")
    SET(ADOLC_WITH_BOOST_ALLOCATOR TRUE)
    SET(_additional_include_dirs OPTIONAL BOOST_INCLUDE_DIRS)
    SET(_additional_library OPTIONAL BOOST_LIBRARIES)
  ELSE()
    SET(ADOLC_WITH_BOOST_ALLOCATOR FALSE)
    SET(_additional_include_dirs)
    SET(_additional_library)
  ENDIF()
ENDIF()


DEAL_II_PACKAGE_HANDLE(ADOLC
  LIBRARIES
    REQUIRED ADOLC_LIBRARY
    ${_additional_library}
  INCLUDE_DIRS 
    REQUIRED ADOLC_INCLUDE_DIR
    ${_additional_include_dirs}
  USER_INCLUDE_DIRS 
    REQUIRED ADOLC_INCLUDE_DIR
    ${_additional_include_dirs}
  CLEAR ADOLC_INCLUDE_DIR ADOLC_LIBRARY ADOLC_SETTINGS_H
    ADOLC_DOUBLE_CAST_CHECK ADOLC_ADOUBLE_OSTREAM_CHECK # clean up checks in configure_adolc.cmake
  )
