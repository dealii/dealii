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

########################################################################
#                                                                      #
#                             Test setup:                              #
#                                                                      #
########################################################################

#
# CTEST_SOURCE_DIRECTORY
# CTEST_BINARY_DIRECTORY
# CTEST_CMAKE_GENERATOR
#
# TRACK
#
# NO_JOBS
# CONFIG_FILE
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
MESSAGE("-- This is CTest ${CTEST_VERSION}")

#
# CTEST_SOURCE_DIR:
#

IF("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_SOURCE_DIRECTORY is not set we just assume that this script
  # was called residing in the source directory.
  #
  SET(CTEST_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
ENDIF()

MESSAGE("-- CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")

#
# CTEST_BINARY_DIR:
#

IF("${CTEST_BINARY_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_BINARY_DIRECTORY is not set we just use the current directory
  # except it is equal to CTEST_SOURCE_DIRECTORY in which case we fail.
  #
  SET(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  IF("${CTEST_BINARY_DIRECTORY}" STREQUAL "${CTEST_SOURCE_DIR}")
    MESSAGE(FATAL_ERROR "
ctest was invoked in the source directory and CTEST_BINARY_DIRECTORY is not set.
Please either call ctest from within a designated build directory, or set CTEST_BINARY_DIRECTORY accordingly."
      )
  ENDIF()
ENDIF()

MESSAGE("-- CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

#
# CTEST_CMAKE_GENERATOR:
#

# Query Generator from build directory (if possible):
IF(EXISTS ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt)
  FILE(STRINGS ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt _generator
    REGEX "^CMAKE_GENERATOR:"
    )
  STRING(REGEX REPLACE "^.*=" "" _generator ${_generator})
ENDIF()

IF("${CTEST_CMAKE_GENERATOR}" STREQUAL "")
  IF(NOT "${_generator}" STREQUAL "")
    SET(CTEST_CMAKE_GENERATOR ${_generator})
  ELSE()
    MESSAGE(FATAL_ERROR "
The build directory is not configured and CTEST_CMAKE_GENERATOR is not set.
different Generator \"${CTEST_CMAKE_GENERATOR}\".
"
     )
  ENDIF()
ELSE()
  # ensure that CTEST_CMAKE_GENERATOR (that was apparantly set) is
  # compatible with the build directory:
  IF( NOT "${CTEST_CMAKE_GENERATOR}" STREQUAL "${_generator}"
      AND NOT "${_generator}" STREQUAL "" )
    MESSAGE(FATAL_ERROR "
The build directory is already set up with Generator \"${_generator}\", but
CTEST_CMAKE_GENERATOR was set to a different Generator \"${CTEST_CMAKE_GENERATOR}\".
"
     )
  ENDIF()
ENDIF()

#
# NO_JOBS:
#

IF("${NO_JOBS}" STREQUAL "")
  SET(NO_JOBS "1")
ENDIF()

MESSAGE("-- NO_JOBS: ${NO_JOBS}")

#
# Finalize:
#

SET(CTEST_NOTES_FILES
  ${CTEST_BINARY_DIRECTORY}/summary.log
  ${CTEST_BINARY_DIRECTORY}/detailed.log
  )


########################################################################
#                                                                      #
#                          Run the testsuite:                          #
#                                                                      #
########################################################################

IF("${TRACK}" STREQUAL "Experimental")

  #
  # TRACK Experimental:
  #
  # It is assumed that the build directory is already configured
  #
  # - Run the configure, build and test stages and submit the results to
  #   the "Experimental" track
  #
  # - No cleanup is done prior to test run
  #

  CTEST_START(Experimental TRACK Experimental)

  IF(NOT "${CONFIG_FILE}" STREQUAL "")
    CTEST_CONFIGURE(OPTIONS -C"${CONFIG_FILE}")
  ELSE()
    CTEST_CONFIGURE()
  ENDIF()

  CTEST_BUILD(TARGET setup_test) # setup tests
  CTEST_BUILD(TARGET) # builds the "all" target

  CTEST_TEST(PARALLEL_LEVEL ${DEAL_II_NO_JOBS})

  CTEST_SUBMIT()

ELSE()

  MESSAGE(FATAL_ERROR "TRACK has to be set TODO")

ENDIF()
