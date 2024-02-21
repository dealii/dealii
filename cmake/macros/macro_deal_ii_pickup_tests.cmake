## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# A macro to set up testing and pick up all tests in the current
# subdirectory.
#
# If TEST_PICKUP_REGEX is set, only tests matching the regex will be
# processed.
#
# Furthermore, the macro sets up (if necessary) deal.II, bash, perl,
# numdiff, and the following variables, that can be overwritten by
# environment or command line:
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
#     deal_ii_pickup_tests()
#

#
# A small helper macro that is used below:
#

macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()


function(pad_string_right output _str _length)
  string(LENGTH "${_str}" _strlen)
  math(EXPR _strlen "${_length} - ${_strlen}")

  if(_strlen GREATER 0)
    if(${CMAKE_VERSION} VERSION_LESS "3.15")
      unset(_pad)
      foreach(_i RANGE 1 ${_strlen}) # inclusive
        string(APPEND _pad " ")
      endforeach()
    else()
      string(REPEAT " " ${_strlen} _pad)
    endif()
    set(_str "${_str}${_pad}")
  endif()

  set(${output} "${_str}" PARENT_SCOPE)
endfunction()


function(pad_string_left output _str _length)
  string(LENGTH "${_str}" _strlen)
  math(EXPR _strlen "${_length} - ${_strlen}")

  if(_strlen GREATER 0)
    if(${CMAKE_VERSION} VERSION_LESS "3.15")
      unset(_pad)
      foreach(_i RANGE 1 ${_strlen}) # inclusive
        string(APPEND _pad " ")
      endforeach()
    else()
      string(REPEAT " " ${_strlen} _pad)
    endif()
    set(_str "${_pad}${_str}")
  endif()

  set(${output} "${_str}" PARENT_SCOPE)
endfunction()



macro(deal_ii_pickup_tests)
  #
  # Initialize two counters to zero:
  #
  set(_number_of_tests 0)
  set(_number_of_test_dependencies 0)

  #
  # Find bash and perl interpreter:
  #

  find_package(UnixCommands REQUIRED)
  find_package(Perl REQUIRED)

  #
  # Ensure that this macro is not called internally:
  #

  if(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    message(FATAL_ERROR
      "\nDEAL_II_PICKUP_TESTS can only be called in external (test sub-) "
      "projects after the inclusion of deal.IIConfig.cmake. It is not "
      "intended for internal use.\n\n"
      )
  endif()

  #
  # Make sure that we know the MPI launcher executable:
  #

  if(${DEAL_II_WITH_MPI})
    if("${DEAL_II_MPIEXEC}" STREQUAL "" OR
       "${DEAL_II_MPIEXEC}" STREQUAL "MPIEXEC_EXECUTABLE-NOTFOUND")
       message(WARNING
         "\nCould not find an MPI launcher program, which is required for "
         "running mpi tests within the testsuite. As a consequence all "
         "tests that require an mpi launcher have been disabled.\n"
         "If you want to run tests with mpi then please configure deal.II "
         "by either setting the MPIEXEC environment variable or the CMake "
         "variable MPIEXEC_EXECUTABLE to a full path to the MPI launcher "
         "program.\n\n"
         )
       set(DEAL_II_WITH_MPI FALSE)
    endif()
  endif()

  #
  # Figure out the category name for the tests:
  #

  if("${ARGN}" STREQUAL "")
    get_filename_component(_category ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  else()
    set(_category "${ARGN}")
  endif()

  #
  # Check that we have the numdiff executable and that it works properly by
  # running a relative tolerance test:
  #

  set_if_empty(NUMDIFF_DIR "$ENV{NUMDIFF_DIR}")
  find_program(NUMDIFF_EXECUTABLE NAMES numdiff
    HINTS ${NUMDIFF_DIR} PATH_SUFFIXES bin
    )
  mark_as_advanced(NUMDIFF_EXECUTABLE)

  set(_file_1 "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/numdiff-1.txt")
  set(_file_2 "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/numdiff-2.txt")
  set(_output "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/numdiff-output.txt")
  file(WRITE "${_file_1}" "0.99999999998\n2.0\n1.0\n")
  file(WRITE "${_file_2}" "1.00000000001\n2.0\n1.0\n")

  set(_command
    "${NUMDIFF_EXECUTABLE}" "-r" "1.0e-8" "--" "${_file_1}" "${_file_2}"
    )

  string(REPLACE ";" " " _contents "${_command}")
  file(WRITE "${_output}" "${_contents}\n")
  execute_process(COMMAND ${_command}
    TIMEOUT 4 # seconds
    OUTPUT_VARIABLE _contents
    ERROR_VARIABLE _contents
    RESULT_VARIABLE _numdiff_tolerance_test_status
    )
  file(APPEND "${_output}" "${_contents}")

  if(NOT "${_numdiff_tolerance_test_status}" STREQUAL "0")
    if(NOT "${_category}" MATCHES "^(quick_tests|performance)$")
      message(WARNING
        "\nCould not find or execute numdiff, which is required for running "
        "most of the tests within the testsuite; failing output has been "
        "recorded in ${_output}"
        "\nAs a consequence all tests that require test output comparison "
        "have been disabled.\nIf you want to run tests that require output "
        "comparison then please specify NUMDIFF_DIR (either as environment "
        "variable or as CMake variable) to a location containing the "
        "binary.\n"
        )
      set(NUMDIFF_EXECUTABLE "")
    endif()
  else()
    file(REMOVE ${_file_1})
    file(REMOVE ${_file_2})
    file(REMOVE ${_output})
  endif()

  #
  # Set various limits:
  #

  set_if_empty(TEST_TIME_LIMIT "$ENV{TEST_TIME_LIMIT}")
  set_if_empty(TEST_TIME_LIMIT 600)

  set_if_empty(TEST_MPI_RANK_LIMIT "$ENV{TEST_MPI_RANK_LIMIT}")
  set_if_empty(TEST_MPI_RANK_LIMIT 0)

  set_if_empty(TEST_THREAD_LIMIT "$ENV{TEST_THREAD_LIMIT}")
  set_if_empty(TEST_THREAD_LIMIT 0)

  #
  # Other variables:
  #

  set_if_empty(TEST_PICKUP_REGEX "$ENV{TEST_PICKUP_REGEX}")

  set_if_empty(ENABLE_PERFORMANCE_TESTS "$ENV{ENABLE_PERFORMANCE_TESTS}")

  set_if_empty(TESTING_ENVIRONMENT "$ENV{TESTING_ENVIRONMENT}")
  set_if_empty(TESTING_ENVIRONMENT "light")

  #
  # ... and finally pick up tests:
  #

  enable_testing()

  set(DEAL_II_SOURCE_DIR) # avoid a bogus warning

  file(GLOB _tests "*.output" "*.run_only")
  foreach(_test ${_tests})
    set(_comparison ${_test})
    get_filename_component(_test ${_test} NAME)

    #
    # Respect TEST_PICKUP_REGEX:
    #

    #
    # Only retain the base name of the test, i.e., remove everything after
    # (and including) the first period ("."):
    #
    string(REGEX REPLACE "\\..*$" "" _regex_name "${_category}/${_test}")
    if( NOT ("${TEST_PICKUP_REGEX}" STREQUAL "" OR
             "${_regex_name}" MATCHES "${TEST_PICKUP_REGEX}" ))
      continue()  # next test
    endif()

    # Disable tests using mpirun if MPI is not enabled
    string(REGEX MATCH "mpirun=" _matches ${_test})
    if (_matches AND NOT DEAL_II_WITH_MPI)
      continue()  # next test
    endif()

    #
    # Query configuration and check whether we support it. Otherwise
    # skip the test.
    #

    set(_op_regex "=|\\.geq\\.|\\.leq\\.|\\.ge\\.|\\.le\\.")

    string(REGEX MATCHALL
      "with_([0-9]|[a-z]|_)*(${_op_regex})(on|off|yes|no|true|false|[0-9]+(\\.[0-9]+)*)"
      _matches ${_test}
      )

    set(_skip_test FALSE)
    foreach(_match ${_matches})
      #
      # Extract feature name, comparison operator, (a possible) boolean and
      # (a possible) version number from the feature constraint:
      #
      string(REGEX REPLACE "^with_(([0-9]|[a-z]|_)*)(${_op_regex}).*" "\\1" _feature ${_match})
      string(TOUPPER ${_feature} _feature)
      string(REGEX MATCH "(${_op_regex})" _operator ${_match})
      string(REGEX REPLACE "^with_(([0-9]|[a-z]|_)*)(${_op_regex}).*$" "\\3" _operator ${_match})
      string(REGEX MATCH "(on|off|yes|no|true|false)$" _boolean ${_match})
      string(REGEX MATCH "([0-9]+(\\.[0-9]+)*)$" _version ${_match})

      #
      # We support two variables: DEAL_II_WITH_<FEATURE> and DEAL_II_<FEATURE>
      #
      set(_variable "DEAL_II_WITH_${_feature}")
      if(NOT DEFINED ${_variable})
        set(_variable "DEAL_II_${_feature}")
        if(NOT DEFINED ${_variable})
          #
          # If a variable is undefined, assume that we cannot configure a
          # given test
          #
          set(_skip_test TRUE)
          continue() # drop out of "foreach(_match ${_matches})"
        endif()
      endif()

      #
      # First process simple yes/no feature constraints:
      #
      if(NOT "${_boolean}" STREQUAL "")
        if(NOT "${_operator}" STREQUAL "=")
          message(FATAL_ERROR "
Invalid syntax in constraint \"${_match}\" in file
\"${_comparison}\":
Comparison operator \"=\" expected for boolean match.\n"
            )
        endif()

        # This is why I hate CMake :-/
        if( (${_variable} AND NOT ${_boolean}) OR
            (NOT ${_variable} AND ${_boolean}) )
          set(_skip_test TRUE)
          continue() # drop out of "foreach(_match ${_matches})"
        endif()
      endif()

      #
      # Process version constraints:
      #
      if(NOT "${_version}" STREQUAL "")

        if( ( NOT ${DEAL_II_WITH_${_feature}} ) OR
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
          set(_skip_test TRUE)
          continue() # drop out of "foreach(_match ${_matches})"
        endif()
      endif()
    endforeach()

    if(_skip_test)
      continue() # next test
    endif()

    #
    # We've made it all the way to here, which means that we actually
    # want to define the test
    #
    string(REGEX REPLACE "\\..*" "" _test ${_test})
    deal_ii_add_test(${_category} ${_test} ${_comparison})

  endforeach()

  if(NOT "${_number_of_tests}" STREQUAL "0")
    # Pad the category to 27 characters -- that's currently the longest
    # category (namely, 'multigrid-global-coarsening'). Pad numbers to 4 digits
    pad_string_right(_category ${_category} 27)
    pad_string_left(_number_of_tests ${_number_of_tests} 4)
    pad_string_left(_number_of_test_dependencies ${_number_of_test_dependencies} 4)
    message(STATUS "Test category ${_category}: ${_number_of_tests} tests (and ${_number_of_test_dependencies} test dependencies)")
  endif()
endmacro()
