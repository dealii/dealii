## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
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

# Note: This file is adapted from https://raw.githubusercontent.com/joaoleal/CppADCodeGen/master/cmake/FindADOLC.cmake

#
# Try to find Adolc
#
# This module exports
#
#   ADOLC_INCLUDE_DIR
#   ADOLC_LIBRARY
#   ADOLC_WITH_ATRIG_ERF
#   ADOLC_WITH_ADVANCED_BRANCHING
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
ENDIF()


DEAL_II_PACKAGE_HANDLE(ADOLC
  LIBRARIES REQUIRED ADOLC_LIBRARY
  INCLUDE_DIRS REQUIRED ADOLC_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED ADOLC_INCLUDE_DIR
  CLEAR ADOLC_INCLUDE_DIR ADOLC_LIBRARY ADOLC_SETTINGS_H
  )
