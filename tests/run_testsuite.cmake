## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


########################################################################
#                                                                      #
#                             Test setup:                              #
#                                                                      #
########################################################################

#
# This is the ctest script for running and submitting build and regression
# tests.
#
# Invoke it in a _build directory_ (or designated build directory) via:
#
#   ctest -S <...>/run_testsuite.cmake
#
# The following configuration variables can be overwritten with
#
#   ctest -D<variable>=<value> [...]
#
#
#   CTEST_SOURCE_DIRECTORY
#     - The source directory of deal.II
#     - If unspecified, "../" relative to the location of this script is
#       used. If this is not a source directory, an error is thrown.
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
#     - If unspecified the current generator of a configured build directory
#       will be used, otherwise "Unix Makefiles".
#
#   TRACK
#     - The track the test should be submitted to. Defaults to
#       "Experimental". Possible values are:
#
#       "Experimental"     - all tests that are not specifically "build" or
#                            "regression" tests should go into this track
#
#       "Build Tests"      - Build tests that configure and build in a clean
#                            directory (without actually running the
#                            testsuite)
#
#       "Regression Tests" - Reserved for "official" regression testers
#
#       "Continuous"       - Reserved for "official" regression testers
#
#   CONFIG_FILE
#     - A configuration file (see ../doc/users/config.sample)
#       that will be used during the configuration stage (invokes
#       # cmake -C ${CONFIG_FILE}). This only has an effect if
#       CTEST_BINARY_DIRECTORY is empty.
#
#   DESCRIPTION
#     - A string that is appended to CTEST_BUILD_NAME
#
#   COVERAGE
#     - If set to ON deal.II will be configured with
#     DEAL_II_SETUP_COVERAGE=ON, CMAKE_BUILD_TYPE=Debug and the
#     ctest_coverage() stage will be run. Test results must go into the
#     "Experimental" section.
#
#   MAKEOPTS
#     - Additional options that will be passed directly to make (or ninja).
#
# Furthermore, the following variables controlling the testsuite can be set
# and will be automatically handed down to cmake:
#
#   NUMDIFF_DIR
#   DIFF_DIR
#   TEST_TIME_LIMIT
#   TEST_PICKUP_REGEX
#
# For details, consult the ./README file.
#

cmake_minimum_required(VERSION 3.15.0)
message("-- This is CTest ${CMAKE_VERSION}")

#
# TRACK: Default to Experimental:
#

if("${TRACK}" STREQUAL "")
  set(TRACK "Experimental")
endif()

if( NOT "${TRACK}" STREQUAL "Experimental"
    AND NOT "${TRACK}" STREQUAL "Build Tests"
    AND NOT "${TRACK}" STREQUAL "Regression Tests"
    AND NOT "${TRACK}" STREQUAL "Continuous" )
  message(FATAL_ERROR "
Unknown TRACK \"${TRACK}\" - see the manual for valid values.
"
    )
endif()

message("-- TRACK:                  ${TRACK}")

#
# CTEST_SOURCE_DIRECTORY:
#

if("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_SOURCE_DIRECTORY is not set we just assume that this script
  # was called residing under ./tests in the source directory
  #
  get_filename_component(CTEST_SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" PATH)

  if(NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CMakeLists.txt)
    message(FATAL_ERROR "
Could not find a suitable source directory. There is no source directory
\"../\" relative to the location of this script. Please, set
CTEST_SOURCE_DIRECTORY manually to the appropriate source directory.
"
      )
  endif()
endif()

message("-- CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")

#
# Read in custom config files:
#

ctest_read_custom_files(${CTEST_SOURCE_DIRECTORY})

#
# CTEST_BINARY_DIRECTORY:
#

if("${CTEST_BINARY_DIRECTORY}" STREQUAL "")
  #
  # If CTEST_BINARY_DIRECTORY is not set we just use the current directory
  # except if it is equal to CTEST_SOURCE_DIRECTORY in which case we fail.
  #
  set(CTEST_BINARY_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  if( ( "${CTEST_BINARY_DIRECTORY}" STREQUAL "${CTEST_SOURCE_DIRECTORY}"
        AND NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/CMakeCache.txt )
      OR "${CTEST_BINARY_DIRECTORY}" STREQUAL "${CMAKE_CURRENT_LIST_DIR}" )
    message(FATAL_ERROR "
ctest was invoked in the source directory or under ./tests and
CTEST_BINARY_DIRECTORY is not set. Please either call ctest from within a
designated build directory, or set CTEST_BINARY_DIRECTORY accordingly.
"
      )
  endif()
endif()

#
# Read in custom config files:
#

ctest_read_custom_files(${CTEST_BINARY_DIRECTORY})

# Make sure that for a build test the directory is empty:
file(GLOB _test ${CTEST_BINARY_DIRECTORY}/*)
if( "${TRACK}" STREQUAL "Build Tests"
    AND NOT "${_test}" STREQUAL "" )
      message(FATAL_ERROR "
TRACK was set to \"Build Tests\" which require an empty build directory.
But files were found in \"${CTEST_BINARY_DIRECTORY}\"
"
        )
endif()

message("-- CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

#
# CTEST_CMAKE_GENERATOR:
#

# Query Generator from build directory (if possible):
if(EXISTS ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt)
  file(STRINGS ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt _generator
    REGEX "^CMAKE_GENERATOR:"
    )
  string(REGEX REPLACE "^.*=" "" _generator ${_generator})
endif()

if("${CTEST_CMAKE_GENERATOR}" STREQUAL "")
  if(NOT "${_generator}" STREQUAL "")
    set(CTEST_CMAKE_GENERATOR ${_generator})
  else()
    # default to "Unix Makefiles"
    set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  endif()
else()
  # ensure that CTEST_CMAKE_GENERATOR (that was apparently set) is
  # compatible with the build directory:
  if( NOT "${CTEST_CMAKE_GENERATOR}" STREQUAL "${_generator}"
      AND NOT "${_generator}" STREQUAL "" )
    message(FATAL_ERROR "
The build directory is already set up with Generator \"${_generator}\", but
CTEST_CMAKE_GENERATOR was set to a different Generator \"${CTEST_CMAKE_GENERATOR}\".
"
     )
  endif()
endif()

message("-- CTEST_CMAKE_GENERATOR:  ${CTEST_CMAKE_GENERATOR}")

#
# CTEST_SITE:
#

if("${CTEST_SITE}" STREQUAL "")
  find_program(HOSTNAME_COMMAND NAMES hostname)
  if(NOT "${HOSTNAME_COMMAND}" MATCHES "-NOTFOUND")
    exec_program(${HOSTNAME_COMMAND} OUTPUT_VARIABLE _hostname)
    string(REGEX REPLACE "\\..*$" "" _hostname ${_hostname})
    set(CTEST_SITE "${_hostname}")
  else()
    # Well, no hostname available. What about:
    set(CTEST_SITE "BobMorane")
  endif()
endif()

message("-- CTEST_SITE:             ${CTEST_SITE}")

if(TRACK MATCHES "^Regression Tests$" AND NOT CTEST_SITE MATCHES "^tester(|-tos|-tng|-ds9)$")
  message(FATAL_ERROR "
I'm sorry ${CTEST_SITE}, I'm afraid I can't do that.
The TRACK \"Regression Tests\" is not for you.
"
    )
endif()

#
# Assemble configuration options, we need it now:
#

if("${MAKEOPTS}" STREQUAL "")
  set(MAKEOPTS $ENV{MAKEOPTS})
endif()

if(NOT "${CONFIG_FILE}" STREQUAL "")
  set(_options "-C${CONFIG_FILE}")
endif()

if("${TRACK}" STREQUAL "Build Tests")
  set(TEST_PICKUP_REGEX "^quick_tests")
endif()

# Pass all relevant variables down to configure:
get_cmake_property(_variables VARIABLES)
foreach(_var ${_variables})
  if( _var MATCHES "^(ENABLE|TEST|TESTING|DEAL_II|ALLOW|WITH|FORCE|COMPONENT)_" OR
      _var MATCHES "^(DOCUMENTATION|EXAMPLES)" OR
      _var MATCHES "^(ADOLC|ARBORX|ARPACK|BOOST|OPENCASCADE|MUPARSER|HDF5|KOKKOS|METIS|MPI)_" OR
      _var MATCHES "^(GINKGO|P4EST|PETSC|SCALAPACK|SLEPC|THREADS|TBB|TRILINOS)_" OR
      _var MATCHES "^(UMFPACK|ZLIB|LAPACK|MUPARSER|MUMPS)_" OR
      _var MATCHES "^(CMAKE|DEAL_II)_(C|CXX|Fortran|BUILD)_(COMPILER|FLAGS)" OR
      _var MATCHES "^CMAKE_BUILD_TYPE$" OR
      _var MATCHES "MAKEOPTS" OR
      ( NOT _var MATCHES "^[_]*CMAKE" AND _var MATCHES "_DIR$" ) )
    list(APPEND _options "-D${_var}=${${_var}}")
  endif()
endforeach()

if(COVERAGE)
  list(APPEND _options "-DDEAL_II_SETUP_COVERAGE=TRUE")
  list(APPEND _options "-DCMAKE_BUILD_TYPE=Debug")
endif()

#
# CTEST_BUILD_NAME:
#

# Append compiler information to CTEST_BUILD_NAME:
if(NOT EXISTS ${CTEST_BINARY_DIRECTORY}/detailed.log)
  # Apparently, ${CTEST_BINARY_DIRECTORY} is not a configured build
  # directory. In this case we need a trick: set up a dummy project and
  # query it for the compiler information.
  file(WRITE ${CTEST_BINARY_DIRECTORY}/query_for_compiler/CMakeLists.txt "
file(WRITE ${CTEST_BINARY_DIRECTORY}/detailed.log
  \"#        CMAKE_CXX_COMPILER:     \${CMAKE_CXX_COMPILER_ID} \${CMAKE_CXX_COMPILER_VERSION} on platform \${CMAKE_SYSTEM_NAME} \${CMAKE_SYSTEM_PROCESSOR}\"
  )"
    )
  execute_process(
    COMMAND ${CMAKE_COMMAND} ${_options} "-G${CTEST_CMAKE_GENERATOR}" .
    OUTPUT_QUIET ERROR_QUIET
    WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}/query_for_compiler
    )
  file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/query_for_compiler)
endif()

if(EXISTS ${CTEST_BINARY_DIRECTORY}/detailed.log)
  file(STRINGS ${CTEST_BINARY_DIRECTORY}/detailed.log _compiler_id
    REGEX "CMAKE_CXX_COMPILER:"
    )
  string(REGEX REPLACE
    "^.*CMAKE_CXX_COMPILER:     \(.*\) on platform.*$" "\\1"
    _compiler_id ${_compiler_id}
    )
  string(REGEX REPLACE "^\(.*\) .*$" "\\1" _compiler_name ${_compiler_id})
  string(REGEX REPLACE "^.* " "" _compiler_version ${_compiler_id})
  string(REGEX REPLACE " " "-" _compiler_id ${_compiler_id})
  if( NOT "${_compiler_id}" STREQUAL "" OR
      _compiler_id MATCHES "CMAKE_CXX_COMPILER" )
    set(CTEST_BUILD_NAME "${_compiler_id}")
  endif()
endif()

#
# Query git information:
#

find_package(Git)

if(NOT GIT_FOUND)
  message(FATAL_ERROR "\nCould not find git. Bailing out.\n"
   )
endif()

execute_process(
   COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:"%h"
   WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
   OUTPUT_VARIABLE _git_WC_INFO
   RESULT_VARIABLE _result
   OUTPUT_STRIP_TRAILING_WHITESPACE
   )

if(NOT ${_result} EQUAL 0)
  message(FATAL_ERROR "\nCould not retrieve git information. Bailing out.\n")
endif()

string(REGEX REPLACE "^\"([^ ]+)\""
         "\\1" _git_WC_SHORTREV "${_git_WC_INFO}")

execute_process(
   COMMAND ${GIT_EXECUTABLE} symbolic-ref HEAD
   WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
   OUTPUT_VARIABLE _git_WC_BRANCH
   RESULT_VARIABLE _result
   OUTPUT_STRIP_TRAILING_WHITESPACE
   )

string(REGEX REPLACE "refs/heads/" ""
  _git_WC_BRANCH "${_git_WC_BRANCH}")

if(NOT "${_git_WC_BRANCH}" STREQUAL "")
  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${_git_WC_BRANCH}")
endif()

#
# Append config file name to CTEST_BUILD_NAME:
#

if(NOT "${CONFIG_FILE}" STREQUAL "")
  get_filename_component(_conf ${CONFIG_FILE} NAME_WE)
  string(REGEX REPLACE "#.*$" "" _conf ${_conf})
  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${_conf}")
endif()

#
# Append DESCRIPTION string to CTEST_BUILD_NAME:
#

if(NOT "${DESCRIPTION}" STREQUAL "")
  set(CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${DESCRIPTION}")
endif()

message("-- CTEST_BUILD_NAME:       ${CTEST_BUILD_NAME}")

#
# Declare files that should be submitted as notes:
#

set(CTEST_NOTES_FILES
  ${CTEST_BINARY_DIRECTORY}/revision.log
  ${CTEST_BINARY_DIRECTORY}/summary.log
  ${CTEST_BINARY_DIRECTORY}/detailed.log
  ${CTEST_BINARY_DIRECTORY}/include/deal.II/base/config.h
  )

#
# Setup coverage:
#

if(COVERAGE)
  if(NOT TRACK MATCHES "Experimental")
    message(FATAL_ERROR "
TRACK must be set to  \"Experimental\" if Coverage is enabled via
COVERAGE=TRUE.
"
      )
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    find_program(GCOV_COMMAND NAMES llvm-cov)
    set(GCOV_COMMAND "${GCOV_COMMAND} gcov")
    if(GCOV_COMMAND MATCHES "-NOTFOUND")
      message(FATAL_ERROR "Coverage enabled but could not find the
                           llvm-cov executable, which is part of LLVM.")
    endif()
  else()
    find_program(GCOV_COMMAND NAMES gcov)
    if(GCOV_COMMAND MATCHES "-NOTFOUND")
      message(FATAL_ERROR "Coverage enabled but could not find the
                           gcov executable, which is part of the
                           GNU Compiler Collection.")
    endif()
  endif()

  set(CTEST_COVERAGE_COMMAND "${GCOV_COMMAND}")
endif()

message("-- COVERAGE:               ${COVERAGE}")


macro(create_targetdirectories_txt)
  #
  # It gets tricky: Fake a TargetDirectories.txt containing _all_ target
  # directories (of the main project and all subprojects) so that the
  # ctest_coverage() actually picks everything up...
  #
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy
    ${CTEST_BINARY_DIRECTORY}/CMakeFiles/TargetDirectories.txt
    ${CTEST_BINARY_DIRECTORY}/CMakeFiles/TargetDirectories.txt.bck
    )
  file(GLOB _subprojects ${CTEST_BINARY_DIRECTORY}/tests/*)
  foreach(_subproject ${_subprojects})
    if(EXISTS ${_subproject}/CMakeFiles/TargetDirectories.txt)
      file(READ ${_subproject}/CMakeFiles/TargetDirectories.txt _var)
      file(APPEND ${CTEST_BINARY_DIRECTORY}/CMakeFiles/TargetDirectories.txt ${_var})
    endif()
  endforeach()
endmacro()

macro(clear_targetdirectories_txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -E rename
    ${CTEST_BINARY_DIRECTORY}/CMakeFiles/TargetDirectories.txt.bck
    ${CTEST_BINARY_DIRECTORY}/CMakeFiles/TargetDirectories.txt
    )
endmacro()

message("-- CMake Options:          ${_options}")

if(NOT "${MAKEOPTS}" STREQUAL "")
  message("-- MAKEOPTS:               ${MAKEOPTS}")
endif()


########################################################################
#                                                                      #
#                          Run the testsuite:                          #
#                                                                      #
########################################################################

# record a status and summary string:
set(_status "neutral")
set(_summary "")

ctest_start(Experimental TRACK ${TRACK})

message("-- Running ctest_update() to query git information")
ctest_update(SOURCE ${CTEST_SOURCE_DIRECTORY})

message("-- Running ctest_configure()")
ctest_configure(OPTIONS "${_options}" RETURN_VALUE _res)

if("${_res}" STREQUAL "0")
  # Only run the build stage if configure was successful:

  message("-- Running ctest_build()")
  set(CTEST_BUILD_FLAGS "${MAKEOPTS}")
  ctest_build(NUMBER_ERRORS _res)

  if("${_res}" STREQUAL "0")
    # Only run tests if the build was successful:

    if(ENABLE_PERFORMANCE_TESTS)
      message("-- Running prune_tests")
      execute_process(COMMAND ${CMAKE_COMMAND}
        --build . --target prune_tests
        -- ${MAKEOPTS}
        WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
        OUTPUT_QUIET
        RESULT_VARIABLE _res
        )

      if(NOT "${_res}" STREQUAL "0")
        message(FATAL_ERROR "
\"prune_tests\" target exited with an error. Bailing out.
"
          )
      endif()

      set(_target setup_tests_performance)
    else()
      set(_target setup_tests)
    endif()

    message("-- Running ${_target}")
    execute_process(COMMAND ${CMAKE_COMMAND}
      --build . --target ${_target}
      -- ${MAKEOPTS}
      WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
      OUTPUT_QUIET
      RESULT_VARIABLE _res
      )

    if(NOT "${_res}" STREQUAL "0")
      message(FATAL_ERROR "
\"setup_tests\" target exited with an error. Bailing out.
"
        )
    endif()

    if(DEAL_II_MSVC)
      set(CTEST_BUILD_CONFIGURATION "${JOB_BUILD_CONFIGURATION}")
    endif()

    message("-- Running ctest_tests()")
    ctest_test()

    if(COVERAGE)
      CREATE_TARGETDIRECTORIES_TXT()
      message("-- Running ctest_coverage()")
      ctest_coverage()
      set (CODE_COV_BASH "${CMAKE_CURRENT_LIST_DIR}/../contrib/utilities/programs/codecov/codecov-bash.sh")
      if (EXISTS ${CODE_COV_BASH})
        message("-- Running codecov-bash")
        execute_process(COMMAND bash "${CODE_COV_BASH}"
                                     "-t ac85e7ce-5316-4bc1-a237-2fe724028c7b" "-x '${GCOV_COMMAND}'"
                        OUTPUT_QUIET)
      endif()
      CLEAR_TARGETDIRECTORIES_TXT()
    endif()

  else()
    # build unsuccessful
    set(_status "failure")
    string(APPEND _summary "\n#   - build failure")
  endif()
else()
  # configure unsuccessful
  set(_status "failure")
  string(APPEND _summary "\n#   - configure failure")
endif()

#
# Create performance test report
#
if(ENABLE_PERFORMANCE_TESTS)
  message("-- Collecting performance measurements\n")
  execute_process(
    COMMAND bash "${CMAKE_CURRENT_LIST_DIR}/performance/collect_measurements" "${CTEST_SITE}"
    WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}
    )
endif()


#
# And finally submit:
#

if(NOT SKIP_SUBMISSION)
  message("-- Running ctest_submit()")
  ctest_submit(RETURN_VALUE _res BUILD_ID _build_id)
  if("${_res}" STREQUAL "0")
    message("-- Submission successful.")
    set(_cdash_url "https://cdash.dealii.org/build/${_build_id}")
  else()
    message("-- Submission failed.")
    set(_cdash_url "-- submission failed --")
  endif()
else()
  set(_cdash_url "-- submission skipped --")
endif()


#
# Grab git revision from our revision.log:
#
file(STRINGS "${CTEST_BINARY_DIRECTORY}/revision.log" _revision REGEX "Revision:")
string(REGEX REPLACE "#.*Revision:  " "" _revision "${_revision}")

#
# Configure or build errors are easy, but determining whether we
# encountered configure or build warnings, or test failures is remarkably
# tricky. None of the ctest_* commands return a value that would help us
# :-(
#
if("${_status}" STREQUAL "neutral")
  #
  # If we made it to this place then configure and build succeeded
  # (otherwise ${_status} would have been set to "failure". So let's try to
  # locate all relevant xml files to query for configure/build warnings and
  # test failures
  #

  # grab tag:
  file(STRINGS ${CTEST_BINARY_DIRECTORY}/Testing/TAG _tag LIMIT_COUNT 1)
  set(_path "${CTEST_BINARY_DIRECTORY}/Testing/${_tag}")
  if(EXISTS "${_path}/Configure.xml" AND EXISTS "${_path}/Build.xml" AND EXISTS "${_path}/Test.xml")
    #
    # All xml files are present. So let's make a decision on "success" or
    # "failure":
    #
    set(_status "success")
    file(STRINGS "${_path}/Configure.xml" _warnings LIMIT_COUNT 1 REGEX "CMake Warning at")
    if(NOT "${_warnings}" STREQUAL "")
      # for the time being configure warnings are not a "failure"
      # condition. This would otherwise create a lot of noise on the
      # regression tester.
      string(APPEND _summary "\n#   - configure warnings")
    endif()
    file(STRINGS "${_path}/Build.xml" _warnings LIMIT_COUNT 1 REGEX "<Warning>")
    if(NOT "${_warnings}" STREQUAL "")
      set(_status "failure")
      string(APPEND _summary "\n#   - build warnings")
    endif()
    file(STRINGS "${_path}/Test.xml" _warnings LIMIT_COUNT 1 REGEX "Status=\"failed\"")
    if(NOT "${_warnings}" STREQUAL "")
      set(_status "failure")
      string(APPEND _summary "\n#   - test failures")
    endif()

    if("${_status}" STREQUAL "success")
      string(APPEND _summary "\n#
# 🎉  🎉  🎉  🎉      Testsuite run succeeded.     🎉  🎉  🎉  🎉")
    endif()
  else()
    message(WARNING "Unable to locate test submission files from TAG.")
  endif()
endif()

message("###
#
# Revision:      ${_revision}
# Site:          ${CTEST_SITE}
# Configuration: ${CTEST_BUILD_NAME}
# CDash URL:     ${_cdash_url}
# Status:        ${_status}${_summary}
#
###")
