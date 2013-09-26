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
# Setup necessary configuration for a testsuite sub project.
#
# A testsuite subproject assumes the following cached variables to be set:
#
#    DEAL_II_BINARY_DIR
#    DEAL_II_SOURCE_DIR
#      - pointing to a source and binary directory of a deal.II build
#
# This file sets up the following options, that can be overwritten by
# environment or command line:
#
#     TEST_DIFF
#     TEST_TIME_LIMIT
#     TEST_PICKUP_REGEX
#

#
# Load all macros:
#

FILE(GLOB _macro_files ${CMAKE_CURRENT_LIST_DIR}/macros/*.cmake)
FOREACH(_file ${_macro_files})
  INCLUDE(${_file})
ENDFOREACH()

#
# Pick up values from environment:
#

SET_IF_EMPTY(DEAL_II_BINARY_DIR $ENV{DEAL_II_BINARY_DIR})
SET_IF_EMPTY(DEAL_II_BINARY_DIR $ENV{DEAL_II_DIR})
SET_IF_EMPTY(DEAL_II_SOURCE_DIR $ENV{DEAL_II_SOURCE_DIR})
SET_IF_EMPTY(TEST_DIFF $ENV{TEST_DIFF})
SET_IF_EMPTY(TEST_TIME_LIMIT $ENV{TEST_TIME_LIMIT})
SET_IF_EMPTY(TEST_PICKUP_REGEX $ENV{TEST_PICKUP_REGEX})

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
  SEPARATE_ARGUMENTS(TEST_DIFF UNIX_COMMAND ${TEST_DIFF})
ENDIF()

#
# Set a default time limit of 600 seconds:
#

SET_IF_EMPTY(TEST_TIME_LIMIT 600)

#
# And finally, enable testing:
#

ENABLE_TESTING()
