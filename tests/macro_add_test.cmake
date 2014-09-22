## ---------------------------------------------------------------------
##
## Copyright (C) 2013, 2014 by the deal.II authors
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
# A Macro to set up tests for the testsuite
#
# The following variables must be set:
#
# TEST_DIFF
#   - specifying the executable and command line of the diff command to use
#
# TEST_TIME_LIMIT
#   - specifying the maximal wall clock time in seconds a test is allowed
#     to run
#
#
# Usage:
#     DEAL_II_ADD_TEST(category test_name comparison_file [ARGN])
#
# This macro assumes that a source file "./tests/category/<test_name>.cc"
# as well as the comparison file "<comparison_file>" is available in the
# testsuite. The output of compiled source file is compared against the
# file comparison file.
#
# This macro gets the following options from the comparison file name (have
# a look at the testsuite documentation for details):
#  - usage of mpirun and number of simultaneous processes
#  - valid build configurations
#  - expected test stage
#

MACRO(DEAL_II_ADD_TEST _category _test_name _comparison_file)

  IF(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    MESSAGE(FATAL_ERROR
      "\nDEAL_II_ADD_TEST can only be called in external test subprojects after "
      "the inclusion of deal.IIConfig.cmake. It is not intended for "
      "internal use.\n\n"
      )
  ENDIF()

  GET_FILENAME_COMPONENT(_file ${_comparison_file} NAME)

  #
  # Determine valid build configurations for this test:
  #

  SET(_configuration)
  IF(_file MATCHES "\\.debug\\.")
    SET(_configuration DEBUG)
  ELSEIF(_file MATCHES "\\.release\\.")
    SET(_configuration RELEASE)
  ENDIF()

  #
  # A "binary" in the output file indicates binary output. In this case we
  # have to switch to plain diff instead of (possibly) numdiff, which can
  # only work on plain text files.
  #
  # TODO: The indirection with ${${_test_diff_variable}} is necessary to
  # avoid quoting issues with the command line :-/ - come up with a fix for
  # that
  #

  SET(_test_diff_variable TEST_DIFF)
  IF(_file MATCHES "\\.binary\\.")
    SET(_test_diff_variable DIFF_EXECUTABLE)
  ENDIF()

  #
  # Determine whether the test should be run with mpirun:
  #

  STRING(REGEX MATCH "mpirun=([0-9]*)" _n_cpu ${_file})
  IF("${_n_cpu}" STREQUAL "")
    SET(_n_cpu 0) # 0 indicates that no mpirun should be used
  ELSE()
    STRING(REGEX REPLACE "^mpirun=([0-9]*)$" "\\1" _n_cpu ${_n_cpu})
  ENDIF()

  #
  # Determine the expected build stage of this test:
  #

  STRING(REGEX MATCH "expect=([a-z]*)" _expect ${_file})
  IF("${_expect}" STREQUAL "")
    SET(_expect "PASSED")
  ELSE()
    STRING(REGEX REPLACE "^expect=([a-z]*)$" "\\1" _expect ${_expect})
    STRING(TOUPPER ${_expect} _expect)
  ENDIF()


  FOREACH(_build ${DEAL_II_BUILD_TYPES})

    ITEM_MATCHES(_match "${_build}" ${_configuration})
    IF(_match OR "${_configuration}" STREQUAL "")

      #
      # Setup a bunch of variables describing the test:
      #
      STRING(TOLOWER ${_build} _build_lowercase)
      SET(_target ${_test_name}.${_build_lowercase}) # target name

      # If _n_cpu is equal to "0", a normal, sequental test will be run,
      # otherwise run the test with mpirun:
      IF("${_n_cpu}" STREQUAL "0")

        SET(_diff_target ${_target}.diff) # diff target name
        SET(_test_full ${_category}/${_test_name}.${_build_lowercase}) # full test name
        SET(_test_directory ${CMAKE_CURRENT_BINARY_DIR}/${_target}) # directory to run the test in
        SET(_run_command ${_target}) # the command to issue

      ELSE()

        SET(_diff_target ${_test_name}.mpirun${_n_cpu}.${_build_lowercase}.diff) # diff target name
        SET(_test_full ${_category}/${_test_name}.mpirun=${_n_cpu}.${_build_lowercase}) # full test name
        SET(_test_directory ${CMAKE_CURRENT_BINARY_DIR}/${_target}/mpirun=${_n_cpu}) # directory to run the test in
        SET(_run_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${_n_cpu} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${_target}/${_target} ${MPIEXEC_POSTFLAGS}) # the command to issue

      ENDIF()

      FILE(MAKE_DIRECTORY ${_test_directory})

      #
      # Add an executable for the current test and set up compile
      # definitions and the full link interface:
      #
      IF(NOT TARGET ${_target})
        # only add the target once

        #
        # Add a "guard file" rule: The purpose of interrupt_guard.cc is to
        # force a complete rerun of this test (BUILD, RUN and DIFF stage)
        # if interrupt_guard.cc is removed by run_test.cmake due to an
        # interruption.
        #
        ADD_CUSTOM_COMMAND(
          OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_target}/interrupt_guard.cc
          COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/${_target}/interrupt_guard.cc
          )

        ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL ${_test_name}.cc
          ${CMAKE_CURRENT_BINARY_DIR}/${_target}/interrupt_guard.cc
          )

        SET_TARGET_PROPERTIES(${_target} PROPERTIES
          LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
          COMPILE_DEFINITIONS "${DEAL_II_USER_DEFINITIONS};${DEAL_II_USER_DEFINITIONS_${_build}}"
          COMPILE_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
          LINKER_LANGUAGE "CXX"
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_target}"
          )
        SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
          INCLUDE_DIRECTORIES "${DEAL_II_INCLUDE_DIRS}"
          )
        SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
          COMPILE_DEFINITIONS
            SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
          )
        TARGET_LINK_LIBRARIES(${_target} ${DEAL_II_TARGET_${_build}})
      ENDIF()

      #
      # Add a top level target to run and compare the test:
      #

      ADD_CUSTOM_COMMAND(OUTPUT ${_test_directory}/output
        COMMAND rm -f ${_test_directory}/failing_output
        COMMAND touch ${_test_directory}/output
        COMMAND
          ${_run_command}
          || (mv ${_test_directory}/output
                 ${_test_directory}/failing_output
              && echo "${_test_full}: BUILD successful."
              && echo "${_test_full}: RUN failed. ------ Result: ${_test_directory}/failing_output"
              && echo "${_test_full}: RUN failed. ------ Partial output:"
              && cat ${_test_directory}/failing_output
              && exit 1)
        COMMAND
          ${PERL_EXECUTABLE} -pi ${DEAL_II_SOURCE_DIR}/tests/normalize.pl
                                 ${_test_directory}/output
        WORKING_DIRECTORY
          ${_test_directory}
        DEPENDS
          ${_target}
          ${DEAL_II_SOURCE_DIR}/tests/normalize.pl
        )
      ADD_CUSTOM_COMMAND(OUTPUT ${_test_directory}/diff
        COMMAND rm -f ${_test_directory}/failing_diff
        COMMAND touch ${_test_directory}/diff
        COMMAND
	  # run diff or numdiff (if available) to determine
	  # whether files are the same. if they are not, output
          # the first few lines of the output of numdiff, followed
          # by the results of regular diff since the latter is just
          # more readable
          ${${_test_diff_variable}} # see comment above about redirection
            ${_comparison_file}
            ${_test_directory}/output
            > ${_test_directory}/diff
          || (mv ${_test_directory}/diff
                 ${_test_directory}/failing_diff
              && echo "${_test_full}: BUILD successful."
              && echo "${_test_full}: RUN successful."
              && echo "${_test_full}: DIFF failed. ------ Source: ${_comparison_file}"
              && echo "${_test_full}: DIFF failed. ------ Result: ${_test_directory}/output"
              && echo "Check ${_test_directory}/output ${_comparison_file}"
              && echo "${_test_full}: DIFF failed. ------ Diff:   ${_test_directory}/failing_diff"
              && echo "${_test_full}: DIFF failed. ------ First 8 lines of numdiff/diff output:"
              && cat ${_test_directory}/failing_diff | head -n 8
              && echo "${_test_full}: DIFF failed. ------ First 50 lines diff output:"
              && ${DIFF_EXECUTABLE} -c ${_comparison_file} ${_test_directory}/output | head -n 50
              && exit 1)
        WORKING_DIRECTORY
          ${_test_directory}
        DEPENDS
          ${_test_directory}/output
          ${_comparison_file}
        )

      ADD_CUSTOM_TARGET(${_diff_target} DEPENDS ${_test_directory}/diff
        COMMAND
             echo "${_test_full}: BUILD successful."
          && echo "${_test_full}: RUN successful."
          && echo "${_test_full}: DIFF successful."
          && echo "${_test_full}: PASSED."
        )

      #
      # And finally add the test:
      #

      ADD_TEST(NAME ${_test_full}
        COMMAND ${CMAKE_COMMAND}
          -DTRGT=${_diff_target}
          -DTEST=${_test_full}
          -DEXPECT=${_expect}
          -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR}
          -DGUARD_FILE=${CMAKE_CURRENT_BINARY_DIR}/${_target}/interrupt_guard.cc
          -P ${DEAL_II_SOURCE_DIR}/tests/run_test.cmake
        WORKING_DIRECTORY ${_test_directory}
        )
      SET_TESTS_PROPERTIES(${_test_full} PROPERTIES
        LABEL "${_category}"
        TIMEOUT ${TEST_TIME_LIMIT}
        )

      #
      # Limit concurrency of mpi tests. We can only set concurrency
      # for the entire test, which includes the compiling and linking
      # stages that are purely sequential. There is no good way to model
      # this without unnecessarily restricting concurrency. Consequently,
      # we just choose to model an "average" concurrency as one half of
      # the number of MPI jobs.
      #
      IF(_n_cpu GREATER 2)
        MATH(EXPR _slots "${_n_cpu} / 2")
        SET_TESTS_PROPERTIES(${_test_full} PROPERTIES PROCESSORS ${_slots})
      ENDIF()

      IF(NOT "${_n_cpu}" STREQUAL "0")
        #
        # We have to be careful not to run different mpirun settings for the
        # same executable in parallel because this triggers a race condition
        # when compiling the not yet existent executable that is shared
        # between the different tests.
        #
        # Luckily CMake has a mechanism to force a test to be run after
        # another has finished (and both are scheduled):
        #
        IF(DEFINED TEST_DEPENDENCIES_${_target})
          SET_TESTS_PROPERTIES(${_test_full} PROPERTIES
            DEPENDS ${TEST_DEPENDENCIES_${_target}}
            )
        ENDIF()
        SET(TEST_DEPENDENCIES_${_target} ${_test_full})
      ENDIF()

    ENDIF()
  ENDFOREACH()
ENDMACRO()
