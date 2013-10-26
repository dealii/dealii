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

#
# Setup necessary configuration in the testsuite subprojects.
# This file is directly included by the test subprojects and not by the
# main project.
#
# It is assumed that the following variables are set:
#
#    DEAL_II_BINARY_DIR
#    DEAL_II_SOURCE_DIR
#      - pointing to a source and binary directory of a deal.II build
#
# This file sets up the following options, that can be overwritten by
# environment or command line:
#
#     TEST_DIFF
#     TEST_OVERRIDE_LOCATION
#     TEST_PICKUP_REGEX
#     TEST_TIME_LIMIT
#

#
# Load all macros:
#
FILE(GLOB _macro_files ${DEAL_II_SOURCE_DIR}/cmake/macros/*.cmake)
FOREACH(_file ${_macro_files})
  INCLUDE(${_file})
ENDFOREACH()

#
# Pick up values from environment:
#
FOREACH(_var
  DEAL_II_BINARY_DIR
  DEAL_II_SOURCE_DIR
  TEST_DIFF
  TEST_TIME_LIMIT
  TEST_PICKUP_REGEX
  TEST_OVERRIDE_LOCATION
  )
  # Environment wins:
  IF(DEFINED ENV{${_var}})
    SET(${_var} $ENV{${_var}})
  ENDIF()
  IF(NOT "${_var}" STREQUAL "")
    SET(${_var} "${${_var}}" CACHE STRING "")
  ENDIF()
ENDFOREACH()

#
# We need deal.II and Perl as external packages:
#
FIND_PACKAGE(deal.II 8.0 REQUIRED
  HINTS ${DEAL_II_BINARY_DIR} ${DEAL_II_DIR}
  )
SET(CMAKE_CXX_COMPILER ${DEAL_II_CXX_COMPILER} CACHE STRING "CXX Compiler.")

FIND_PACKAGE(Perl REQUIRED)

#
# We need a diff tool, preferably numdiff:
#
FIND_PROGRAM(DIFF_EXECUTABLE
  NAMES diff
  )

FIND_PROGRAM(NUMDIFF_EXECUTABLE
  NAMES numdiff
  HINTS ${NUMDIFF_DIR}
  PATH_SUFFIXES bin
  )

MARK_AS_ADVANCED(DIFF_EXECUTABLE NUMDIFF_EXECUTABLE)

IF( NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND"
    AND DIFF_EXECUTABLE MATCHES "-NOTFOUND" )
  MESSAGE(FATAL_ERROR
    "Could not find diff or numdiff. One of those are required for running the testsuite."
    )
ENDIF()

IF("${TEST_DIFF}" STREQUAL "")
  IF(NOT NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND")
      SET(TEST_DIFF ${NUMDIFF_EXECUTABLE} -a 1e-6 -s ' \\t\\n:')
  ELSE()
      SET(TEST_DIFF ${DIFF_EXECUTABLE})
  ENDIF()
ELSE()
  # TODO: I have no idea how to prepare a custom string comming possibly
  # through two layers of command line into a list...
  SEPARATE_ARGUMENTS(TEST_DIFF ${TEST_DIFF})
ENDIF()

#
# Set a default time limit of 600 seconds:
#
SET_IF_EMPTY(TEST_TIME_LIMIT 600)

#
# And finally, enable testing:
#
ENABLE_TESTING()

#
# A custom target that does absolutely nothing. It is used in the main
# project to trigger a "make rebuild_cache" if necessary.
#
ADD_CUSTOM_TARGET(regenerate)
