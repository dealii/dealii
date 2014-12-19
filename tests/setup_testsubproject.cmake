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
#     TEST_PICKUP_REGEX
#     TEST_TIME_LIMIT
#

#
# Silence warnings:
#
FOREACH(_var
    MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_POSTFLAGS MPIEXEC_PREFLAGS
    )
  SET(${_var} ${${_var}})
ENDFOREACH()


#
# Load all macros:
#
FILE(GLOB _macro_files ${DEAL_II_SOURCE_DIR}/cmake/macros/*.cmake)
FOREACH(_file ${_macro_files})
  INCLUDE(${_file})
ENDFOREACH()

INCLUDE(${DEAL_II_SOURCE_DIR}/tests/macro_add_test.cmake)
INCLUDE(${DEAL_II_SOURCE_DIR}/tests/macro_pickup_tests.cmake)

#
# Pick up values from environment:
#
FOREACH(_var
  DEAL_II_BINARY_DIR
  DEAL_II_SOURCE_DIR
  TEST_DIFF
  TEST_TIME_LIMIT
  TEST_PICKUP_REGEX
  )
  # Environment wins:
  IF(DEFINED ENV{${_var}})
    SET(${_var} $ENV{${_var}})
  ENDIF()
  IF(NOT "${_var}" STREQUAL "")
    SET(${_var} "${${_var}}" CACHE STRING "" FORCE)
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
  HINTS ${DIFF_DIR}
  PATH_SUFFIXES bin
  )

FIND_PROGRAM(NUMDIFF_EXECUTABLE
  NAMES numdiff
  HINTS ${NUMDIFF_DIR}
  PATH_SUFFIXES bin
  )

MARK_AS_ADVANCED(DIFF_EXECUTABLE NUMDIFF_EXECUTABLE)

IF("${TEST_DIFF}" STREQUAL "")
  #
  # No TEST_DIFF is set, specify one:
  #

  IF(NOT NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND")
    SET(TEST_DIFF ${NUMDIFF_EXECUTABLE} -a 1e-6 -r 1e-8 -s ' \\t\\n:<>=,;')
    IF(DIFF_EXECUTABLE MATCHES "-NOTFOUND")
      SET(DIFF_EXECUTABLE ${NUMDIFF_EXECUTABLE})
    ENDIF()
  ELSEIF(NOT DIFF_EXECUTABLE MATCHES "-NOTFOUND")
    SET(TEST_DIFF ${DIFF_EXECUTABLE})
  ELSE()
    MESSAGE(FATAL_ERROR
      "Could not find diff or numdiff. One of those are required for running the testsuite.\n"
      "Please specify TEST_DIFF by hand."
      )
  ENDIF()
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
