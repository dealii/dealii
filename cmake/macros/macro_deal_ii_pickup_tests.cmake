## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2022 by the deal.II authors
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

#
# A macro to set up testing and pick up all tests in the current
# subdirectory.
#
# If TEST_PICKUP_REGEX is set, only tests matching the regex will be
# processed.
#
# Furthermore, the macro sets up (if necessary) deal.II, perl, numdiff,
# and the following variables, that can be overwritten by environment or
# command line:
#
#     TEST_LIBRARIES
#     TEST_LIBRARIES_DEBUG
#     TEST_LIBRARIES_RELEASE
#       - Specify additional libraries (and targets) to link against.
#
#     TEST_TARGET or
#     TEST_TARGET_DEBUG and TEST_TARGET_RELEASE
#       - Specifies a test target to be executed for a parameter run.
#
#     TEST_TIME_LIMIT
#       - Specifies the maximal wall clock time in seconds a test is
#         allowed to run. Defaults to 600.
#     TEST_MPI_RANK_LIMIT
#       - Specifies the maximal number of MPI ranks that can be used. If a
#         test variant configures a larger number of MPI ranks (via
#         .mpirun=N. in the output file) than this limit the test will be
#         dropped. The special value 0 enforces no limit. Defaults to 0.
#     TEST_THREAD_LIMIT
#       - Specifies the maximal number of worker threads that can should be
#         used by the threading backend. If a test variant configures a
#         larger number of threads (via .threads=N. in the output file)
#         than this limit the test will be dropped. Note that individual
#         tests might exceed this limit by calling
#         MultithreadInfo::set_thread_limit(), or by manually creating
#         additional threads. The special value 0 enforces no limit.
#         Defaults to 0.
#
#     TEST_PICKUP_REGEX
#       - A regular expression to select only a subset of tests during setup.
#         An empty string is interpreted as a catchall (this is the default).
#
#     ENABLE_PERFORMANCE_TESTS
#       - If defined and set to true the execution of performance tests
#         will be enabled.
#
#     TESTING_ENVIRONMENT
#       - Specifies the performance test testing environment. Valid options
#         are:
#          * "light":  mobile laptop, >=2 physical cores, >=8GB RAM
#          * "medium": workstation, >=8 physical cores, >=32GB RAM
#          * "heavy":  compute node, >=32 physical cores, >=128GB RAM
#
# numdiff is used for the comparison of test results. Its location can be
# specified with NUMDIFF_DIR.
#
# Usage:
#     DEAL_II_PICKUP_TESTS()
#

#
# Two very small macros that are used below:
#

MACRO(SET_IF_EMPTY _variable)
  IF("${${_variable}}" STREQUAL "")
    SET(${_variable} ${ARGN})
  ENDIF()
ENDMACRO()

MACRO(ITEM_MATCHES _var _regex)
  SET(${_var})
  FOREACH (_item ${ARGN})
    IF("${_item}" MATCHES ${_regex})
      SET(${_var} TRUE)
      BREAK()
    ENDIF()
  ENDFOREACH()
ENDMACRO()


MACRO(DEAL_II_PICKUP_TESTS)

  IF(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    MESSAGE(FATAL_ERROR
      "\nDEAL_II_PICKUP_TESTS can only be called in external (test sub-) "
      "projects after the inclusion of deal.IIConfig.cmake. It is not "
      "intended for internal use.\n\n"
      )
  ENDIF()

  #
  # Necessary external interpreters and programs:
  #
  IF(${DEAL_II_WITH_MPI})
    IF("${DEAL_II_MPIEXEC}" STREQUAL "" OR
       "${DEAL_II_MPIEXEC}" STREQUAL "MPIEXEC_EXECUTABLE-NOTFOUND")
      MESSAGE(FATAL_ERROR "Could not find an MPI launcher program, which is required "
"for running the testsuite. Please explicitly specify MPIEXEC_EXECUTABLE to CMake "
"as a full path to the MPI launcher program.")
    ENDIF()
  ENDIF()

  IF(DEAL_II_WITH_CUDA)
    FIND_PACKAGE(CUDA)
    SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -std=c++14 -arch=sm_35 -Xcompiler ${OpenMP_CXX_FLAGS})
  ENDIF()

  FIND_PACKAGE(Perl REQUIRED)

  FIND_PROGRAM(NUMDIFF_EXECUTABLE
    NAMES numdiff
    HINTS ${NUMDIFF_DIR}
    PATH_SUFFIXES bin
    )

  MARK_AS_ADVANCED(NUMDIFF_EXECUTABLE)

  IF( "${NUMDIFF_EXECUTABLE}" MATCHES "NUMDIFF_EXECUTABLE-NOTFOUND")
    MESSAGE(FATAL_ERROR
      "Could not find numdiff, which is required for running the testsuite.\n"
      "Please specify NUMDIFF_DIR to a location containing the binary."
      )
  ENDIF()

  #
  # Check that numdiff can run and terminate successfully:
  #
  EXECUTE_PROCESS(COMMAND ${NUMDIFF_EXECUTABLE} "-v"
    TIMEOUT 4 # seconds
    OUTPUT_QUIET
    ERROR_QUIET
    RESULT_VARIABLE _diff_program_status
    )

  IF(NOT "${_diff_program_status}" STREQUAL "0")
    MESSAGE(FATAL_ERROR
      "\nThe command \"${NUMDIFF_EXECUTABLE} -v\" did not run correctly: it "
      "either failed to exit after a few seconds or returned a nonzero exit "
      "code. The test suite cannot be set up without this program, so please "
      "reinstall it and then run the test suite setup command again.\n")
  ENDIF()

  #
  # Also check that numdiff is not a symlink to diff by running a relative
  # tolerance test.
  #
  STRING(FIND "${NUMDIFF_EXECUTABLE}" "numdiff" _found_numdiff_binary)
  IF(NOT "${_found_numdiff_binary}" STREQUAL "-1")
    STRING(RANDOM _suffix)
    SET(_first_test_file_name
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/numdiff-test-${_suffix}-1.txt")
    SET(_second_test_file_name
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/numdiff-test-${_suffix}-2.txt")
    FILE(WRITE "${_first_test_file_name}" "0.99999999998\n2.0\n1.0\n")
    FILE(WRITE "${_second_test_file_name}" "1.00000000001\n2.0\n1.0\n")

    EXECUTE_PROCESS(COMMAND ${NUMDIFF_EXECUTABLE}
      "-r" "1.0e-8" "--" "${_first_test_file_name}" "${_second_test_file_name}"
      TIMEOUT 4 # seconds
      OUTPUT_QUIET
      ERROR_QUIET
      RESULT_VARIABLE _numdiff_tolerance_test_status
      )

    #
    # Tidy up:
    #
    FILE(REMOVE ${_first_test_file_name})
    FILE(REMOVE ${_second_test_file_name})

    IF(NOT "${_numdiff_tolerance_test_status}" STREQUAL "0")
      MESSAGE(FATAL_ERROR
        "\nThe detected numdiff executable was not able to pass a simple "
        "relative tolerance test. This usually means that either numdiff "
        "was misconfigured or that it is a symbolic link to diff. "
        "The test suite needs numdiff to work correctly: please reinstall "
        "numdiff and run the test suite configuration again.\n")
    ENDIF()
  ENDIF()

  #
  # Set various limits:
  #

  SET_IF_EMPTY(TEST_TIME_LIMIT "$ENV{TEST_TIME_LIMIT}")
  SET_IF_EMPTY(TEST_TIME_LIMIT 600)

  SET_IF_EMPTY(TEST_MPI_RANK_LIMIT "$ENV{TEST_MPI_RANK_LIMIT}")
  SET_IF_EMPTY(TEST_MPI_RANK_LIMIT 0)

  SET_IF_EMPTY(TEST_THREAD_LIMIT "$ENV{TEST_THREAD_LIMIT}")
  SET_IF_EMPTY(TEST_THREAD_LIMIT 0)

  #
  # Other variables:
  #

  SET_IF_EMPTY(TEST_PICKUP_REGEX "$ENV{TEST_PICKUP_REGEX}")

  SET_IF_EMPTY(ENABLE_PERFORMANCE_TESTS "$ENV{ENABLE_PERFORMANCE_TESTS}")

  SET_IF_EMPTY(TESTING_ENVIRONMENT "$ENV{TESTING_ENVIRONMENT}")
  SET_IF_EMPTY(TESTING_ENVIRONMENT "light")

  #
  # ... and finally pick up tests:
  #

  ENABLE_TESTING()

  IF("${ARGN}" STREQUAL "")
    GET_FILENAME_COMPONENT(_category ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  ELSE()
    SET(_category "${ARGN}")
  ENDIF()

  SET(DEAL_II_SOURCE_DIR) # avoid a bogus warning

  FILE(GLOB _tests "*.output" "*.run_only")
  FOREACH(_test ${_tests})
    SET(_comparison ${_test})
    GET_FILENAME_COMPONENT(_test ${_test} NAME)

    #
    # Respect TEST_PICKUP_REGEX:
    #

    #
    # Only retain the base name of the test, i.e., remove everything after
    # (and including) the first period ("."):
    #
    STRING(REGEX REPLACE "\\..*$" "" _regex_name "${_category}/${_test}")
    IF( NOT ("${TEST_PICKUP_REGEX}" STREQUAL "" OR
             "${_regex_name}" MATCHES "${TEST_PICKUP_REGEX}" ))
      CONTINUE()  # next test
    ENDIF()

    # Disable tests using mpirun if MPI is not enabled
    STRING(REGEX MATCH "mpirun=" _matches ${_test})
    IF (_matches AND NOT DEAL_II_WITH_MPI)
      CONTINUE()  # next test
    ENDIF()

    #
    # Query configuration and check whether we support it. Otherwise
    # skip the test.
    #

    SET(_op_regex "=|\\.geq\\.|\\.leq\\.|\\.ge\\.|\\.le\\.")

    STRING(REGEX MATCHALL
      "with_([0-9]|[a-z]|_)*(${_op_regex})(on|off|yes|no|true|false|[0-9]+(\\.[0-9]+)*)"
      _matches ${_test}
      )

    SET(_skip_test FALSE)
    FOREACH(_match ${_matches})
      #
      # Extract feature name, comparison operator, (a possible) boolean and
      # (a possible) version number from the feature constraint:
      #
      STRING(REGEX REPLACE "^with_(([0-9]|[a-z]|_)*)(${_op_regex}).*" "\\1" _feature ${_match})
      STRING(TOUPPER ${_feature} _feature)
      STRING(REGEX MATCH "(${_op_regex})" _operator ${_match})
      STRING(REGEX REPLACE "^with_(([0-9]|[a-z]|_)*)(${_op_regex}).*$" "\\3" _operator ${_match})
      STRING(REGEX MATCH "(on|off|yes|no|true|false)$" _boolean ${_match})
      STRING(REGEX MATCH "([0-9]+(\\.[0-9]+)*)$" _version ${_match})

      #
      # We support two variables: DEAL_II_WITH_<FEATURE> and DEAL_II_<FEATURE>
      #
      SET(_variable "DEAL_II_WITH_${_feature}")
      IF(NOT DEFINED ${_variable})
        SET(_variable "DEAL_II_${_feature}")
        IF(NOT DEFINED ${_variable})
          #
          # If a variable is undefined, assume that we cannot configure a
          # given test
          #
          SET(_skip_test TRUE)
          CONTINUE() # drop out of "FOREACH(_match ${_matches})"
        ENDIF()
      ENDIF()

      #
      # First process simple yes/no feature constraints:
      #
      IF(NOT "${_boolean}" STREQUAL "")
        IF(NOT "${_operator}" STREQUAL "=")
          MESSAGE(FATAL_ERROR "
Invalid syntax in constraint \"${_match}\" in file
\"${_comparison}\":
Comparison operator \"=\" expected for boolean match.\n"
            )
        ENDIF()

        # This is why I hate CMake :-/
        IF( (${_variable} AND NOT ${_boolean}) OR
            (NOT ${_variable} AND ${_boolean}) )
          SET(_skip_test TRUE)
          CONTINUE() # drop out of "FOREACH(_match ${_matches})"
        ENDIF()
      ENDIF()

      #
      # Process version constraints:
      #
      IF(NOT "${_version}" STREQUAL "")

        IF( ( NOT ${DEAL_II_WITH_${_feature}} ) OR
            ( "${_operator}" STREQUAL "=" AND
              NOT "${DEAL_II_${_feature}_VERSION}" VERSION_EQUAL "${_version}" ) OR
            ( "${_operator}" STREQUAL ".ge." AND
              NOT "${DEAL_II_${_feature}_VERSION}" VERSION_GREATER "${_version}" ) OR
            ( "${_operator}" STREQUAL ".le." AND
              NOT "${DEAL_II_${_feature}_VERSION}" VERSION_LESS "${_version}" ) OR
            ( "${_operator}" STREQUAL ".geq." AND
              "${DEAL_II_${_feature}_VERSION}" VERSION_LESS "${_version}" ) OR
            ( "${_operator}" STREQUAL ".leq." AND
              "${DEAL_II_${_feature}_VERSION}" VERSION_GREATER "${_version}" ) )
          SET(_skip_test TRUE)
          CONTINUE() # drop out of "FOREACH(_match ${_matches})"
        ENDIF()
      ENDIF()
    ENDFOREACH()

    IF(_skip_test)
      CONTINUE() # next test
    ENDIF()

    #
    # We've made it all the way to here, which means that we actually
    # want to define the test
    #
    STRING(REGEX REPLACE "\\..*" "" _test ${_test})
    DEAL_II_ADD_TEST(${_category} ${_test} ${_comparison})

  ENDFOREACH()
ENDMACRO()
