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
# This is the test script for running the testsuite.
#
# Invoke it in a _build directory_ (or designated build directory) via:
#
#   ctest -S <...>/run_testsuite.cmake
#
# The following configuration variables can be overriden with
#
#   ctest -D<variable>=<value> [...]
#
#
#   CTEST_SOURCE_DIRECTORY
#     - The source directory of deal.II (usually ending in "[...]/deal.II"
#       (equivalent to https://svn.dealii.org/trunk/deal.II)
#       Note: This is _not_ the test directory ending in "[...]/tests"
#     - If unspecified, "../deal.II" relative to the location of this
#       script is used. If this is not a source directory, an error is
#       thrown.
#
#   CTEST_BINARY_DIRECTORY
#     - The designated build directory (already configured, empty, or non
#       existent - see the information about TRACKs what will happen)
#     - If unspecified the current directory is used. If the current
#       directory is equal to CTEST_SOURCE_DIRECTORY or the "tests"
#       directory, an error is thrown.
#
#   CTEST_CMAKE_GENERATOR
#     - The CMake Generator to use (e.g. "Unix Makefiles", or "Ninja", see
#       $ man cmake)
#     - If unspecified the generator of a configured build directory will
#       be used, otherwise "Unix Makefiles".
#
#   TRACK
#     - TODO (defaults to "Experimental")
#
#   CONFIG_FILE
#     - A configuration file (see ../deal.II/docs/development/Config.sample)
#       that will be used during the configuration stage (invokes
#       # cmake -C ${CONFIG_FILE}). This only has an effect if
#       CTEST_BINARY_DIRECTORY is empty.
#
# Furthermore, the following variables controlling the testsuite can be set
# and will be automatically handed down to cmake:
#
#   TEST_DIFF
#   TEST_TIME_LIMIT
#   TEST_PICKUP_REGEX
#   NUMDIFF_DIR
#
# For details, consult the ./README file.
#


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
MESSAGE("-- This is CTest ${CTEST_VERSION}")

#
# TRACK: Default to Experimental:
#

IF("${TRACK}" STREQUAL "")
  SET(TRACK "Experimental")
ENDIF()

IF( NOT "${TRACK}" STREQUAL "Experimental"
    AND NOT "${TRACK}" STREQUAL "Build Tests" )
  MESSAGE(FATAL_ERROR "TODO: Unknown TRACK \"${TRACK}\"")
ENDIF()

MESSAGE("-- TRACK:                  ${TRACK}")

#
# CTEST_SOURCE_DIRECTORY:
#

IF("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_SOURCE_DIRECTORY is not set we just assume that this script
  # was called residing under ../tests relative to the source directory.
  #
  GET_FILENAME_COMPONENT(_path "${CMAKE_CURRENT_LIST_DIR}" PATH)
  SET(CTEST_SOURCE_DIRECTORY ${_path}/deal.II)

  IF(NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CMakeLists.txt)
    MESSAGE(FATAL_ERROR "
Could not find a suitable source directory. Either the run_testsuite.cmake
script was called residing outside of \"tests\", or there is no source
directory \"../deal.II\" relative to the tests directory.
Please, set CTEST_SOURCE_DIRECTORY manually to the appropriate source
directory.
"
      )
  ENDIF()
ENDIF()

MESSAGE("-- CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")

#
# CTEST_BINARY_DIRECTORY:
#

IF("${CTEST_BINARY_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_BINARY_DIRECTORY is not set we just use the current directory
  # except if it is equal to CTEST_SOURCE_DIRECTORY in which case we fail.
  #
  SET(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  IF( "${CTEST_BINARY_DIRECTORY}" STREQUAL "${CTEST_SOURCE_DIR}"
      OR "${CTEST_BINARY_DIRECTORY}" STREQUAL "${CMAKE_CURRENT_LIST_DIR}" )
    MESSAGE(FATAL_ERROR "
ctest was invoked in the source directory (or test source directory) and CTEST_BINARY_DIRECTORY is not set.
Please either call ctest from within a designated build directory, or set CTEST_BINARY_DIRECTORY accordingly.
"
      )
  ENDIF()
ENDIF()

# Make sure that for a build test the directory is empty:
FILE(GLOB _test ${CTEST_BINARY_DIRECTORY}/*)
IF( "${TRACK}" STREQUAL "Build Tests"
    AND NOT "${_test}" STREQUAL "" )
    MESSAGE(FATAL_ERROR "
TRACK was set to \"Build Tests\" which require an empty build directory.
But files were found in \"${CTEST_BINARY_DIRECTORY}\"
"
      )
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
    # default to "Unix Makefiles"
    SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")
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
# CTEST_UPDATE_COMMAND:
#

FIND_PACKAGE(Subversion QUIET)
IF(SUBVERSION_FOUND)
  SET(CTEST_UPDATE_COMMAND ${Subversion_SVN_EXECUTABLE})
ENDIF()

MESSAGE("-- CTEST_UPDATE_COMMAND:   ${CTEST_UPDATE_COMMAND}")

#
# CTEST_SITE:
#

FIND_PROGRAM(HOSTNAME_COMMAND NAMES hostname)
EXEC_PROGRAM(${HOSTNAME_COMMAND} OUTPUT_VARIABLE _hostname)
SET(CTEST_SITE "${_hostname}")

MESSAGE("-- CTEST_SITE:             ${CTEST_SITE}")

#
# Assemble configuration options, we need it now:
#

IF(NOT "${CONFIG_FILE}" STREQUAL "")
  SET(_options "-C${CONFIG_FILE}")
ENDIF()

IF("${TRACK}" STREQUAL "Build Tests")
  SET(TEST_PICKUP_REGEX "^build_tests")
ENDIF()

# Pass all relevant "TEST_" variables down to configure:
GET_CMAKE_PROPERTY(_variables VARIABLES)
FOREACH(_var ${_variables})
  IF(_var MATCHES
      "^(TEST_DIFF|TEST_TIME_LIMIT|TEST_PICKUP_REGEX|NUMDIFF_DIR)$"
      )
    LIST(APPEND _options "-D${_var}=${${_var}}")
  ENDIF()
ENDFOREACH()

#
# CTEST_BUILD_NAME:
#

SET(CTEST_BUILD_NAME
  "${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM_NAME}"
  )

#
# Append compiler information to CTEST_BUILD_NAME:
#

IF(NOT EXISTS ${CTEST_BINARY_DIRECTORY}/detailed.log)
  # Apparently, ${CTEST_BINARY_DIRECTORY} is not a configured build
  # directory. In this case we need a trick: set up a dummy project and
  # query it for the compiler information.
  FILE(WRITE ${CTEST_BINARY_DIRECTORY}/query_for_compiler/CMakeLists.txt "
FILE(WRITE ${CTEST_BINARY_DIRECTORY}/detailed.log
  \"#        CMAKE_CXX_COMPILER:     \${CMAKE_CXX_COMPILER_ID} \${CMAKE_CXX_COMPILER_VERSION} on platform \${CMAKE_SYSTEM_NAME} \${CMAKE_SYSTEM_PROCESSOR}\"
  )"
    )
  EXECUTE_PROCESS(
    COMMAND ${CMAKE_COMMAND} ${_options} "-G${CTEST_CMAKE_GENERATOR}" .
    OUTPUT_QUIET ERROR_QUIET
    WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}/query_for_compiler
    )
  FILE(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/query_for_compiler)
ENDIF()

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

#
# Append subversion information to CTEST_BUILD_NAME:
#

IF(NOT TRACK MATCHES "Foobar" AND SUBVERSION_FOUND) #TODO: Exclude Continuous TRACKS
  EXECUTE_PROCESS(
    COMMAND ${Subversion_SVN_EXECUTABLE} info ${CTEST_SOURCE_DIRECTORY}
    OUTPUT_QUIET ERROR_QUIET
    RESULT_VARIABLE _result
    )
  IF(${_result} EQUAL 0)
    Subversion_WC_INFO(${CTEST_SOURCE_DIRECTORY} _svn)
    STRING(REGEX REPLACE "^${_svn_WC_ROOT}/" "" _branch ${_svn_WC_URL})
    STRING(REGEX REPLACE "^branches/" "" _branch ${_branch})
    STRING(REGEX REPLACE "/deal.II$" "" _branch ${_branch})
    SET(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${_branch}-r${_svn_WC_REVISION}")
  ENDIF()
ENDIF()

MESSAGE("-- CTEST_BUILD_NAME:       ${CTEST_BUILD_NAME}")

#
# Write revision log:
#

IF(DEFINED _svn_WC_REVISION)
  FILE(WRITE ${CTEST_BINARY_DIRECTORY}/revision.log
"###
#
#  SVN information:
#        SVN_WC_URL:               ${_svn_WC_URL}
#        SVN_WC_REVISION:          ${_svn_WC_REVISION}
#        SVN_WC_LAST_CHANGED_DATE: ${_svn_WC_LAST_CHANGED_DATE}
#
###"
    )
ELSE()
  FILE(WRITE ${CTEST_BINARY_DIRECTORY}/revision.log
"###
#
#  No SVN information available.
#
###"
    )
ENDIF()

#
# Declare files that should be submitted as notes:
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

CTEST_START(Experimental TRACK ${TRACK})

CTEST_CONFIGURE(OPTIONS "${_options}")

CTEST_BUILD(TARGET) # run all target

# TODO: Run this during the BUILD stage...
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
  --build ${CTEST_BINARY_DIRECTORY} --target setup_test
  )

CTEST_TEST()

CTEST_SUBMIT()
