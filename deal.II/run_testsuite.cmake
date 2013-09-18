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
# Configuration options:
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
# CTEST_SITE:
#

FIND_PROGRAM(HOSTNAME_COMMAND NAMES hostname)
EXEC_PROGRAM(${HOSTNAME_COMMAND} OUTPUT_VARIABLE _hostname)
SET(CTEST_SITE "${_hostname}")

#
# CTEST_UPDATE_COMMAND
#

FIND_PACKAGE(Subversion QUIET)
IF(SUBVERSION_FOUND)
  SET(CTEST_UPDATE_COMMAND ${Subversion_SVN_EXECUTABLE})
ENDIF()

#
# CTEST_BUILD_NAME
#

SET(CTEST_BUILD_NAME
  "${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}"
  )

# Append compiler information if available...
IF(EXISTS ${CTEST_BINARY_DIRECTORY}/detailed.log)
  FILE(STRINGS ${CTEST_BINARY_DIRECTORY}/detailed.log _compiler_id
    REGEX "CMAKE_CXX_COMPILER:"
    )
  STRING(REGEX REPLACE
    "^.*CMAKE_CXX_COMPILER:     \(.*\) on platform.*$" "\\1"
    _compiler_id ${_compiler_id}
    )
  STRING(REGEX REPLACE " " "-" _compiler_id ${_compiler_id})
  IF( NOT "${_compiler_id}" STREQUAL "" OR
      _compiler_id MATCHES "CMAKE_CXX_COMPILER" )
    SET(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${_compiler_id}")
  ENDIF()
ENDIF()

IF(NOT TRACK MATCHES "Foobar") #TODO: Exclude Continuous TRACKS
  #Append subversion revision:
  IF(SUBVERSION_FOUND)
    EXECUTE_PROCESS(COMMAND ${Subversion_SVN_EXECUTABLE} info
      ${CTEST_SOURCE_DIRECTORY} OUTPUT_VARIABLE _tmp RESULT_VARIABLE _result
      )
    IF(${_result} EQUAL 0)
      Subversion_WC_INFO(${CTEST_SOURCE_DIRECTORY} _svn)
      SET(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-r${_svn_WC_REVISION}")
    ENDIF()
  ENDIF()
ENDIF()

#
# Finalize:
#

SET(CTEST_NOTES_FILES
  ${CTEST_BINARY_DIRECTORY}/revision.log
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
  # - No cleanup is done
  #

  CTEST_START(Experimental TRACK Experimental)

  IF(NOT "${CONFIG_FILE}" STREQUAL "")
    CTEST_CONFIGURE(OPTIONS -C"${CONFIG_FILE}")
  ELSE()
    CTEST_CONFIGURE()
  ENDIF()

  CTEST_BUILD(TARGET) # run all target

  # TODO: Run this during the BUILD stage...
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
    --build ${CTEST_BINARY_DIRECTORY} --target setup_test
    )

  CTEST_TEST(PARALLEL_LEVEL ${NO_JOBS})

  CTEST_SUBMIT()

ELSE()

  MESSAGE(FATAL_ERROR "TRACK has to be set TODO")

ENDIF()
