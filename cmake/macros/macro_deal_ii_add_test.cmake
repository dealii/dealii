## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2019 by the deal.II authors
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
# A Macro to set up tests for the testsuite
#
#
# The testsuite distinguishes two different kinds of tests:
#
# - A combination of a source file "${test_name}.cc" (containing a main
#   function) with a file "${comparison_file}" defines an executable that
#   is compiled and linked against deal.II. Its output is compared with the
#   comparison file. Additional libraries (like a library from a user
#   project with code to test) the target should be linked against can be
#   specified by
#
#     TEST_LIBRARIES
#     TEST_LIBRARIES_DEBUG
#     TEST_LIBRARIES_RELEASE
#
# - A combination of a parameter file "${test_name}.prm" with a file
#   "${comparison_file}" describes the configuration of an already compiled
#   executable that should just be run with its output being compared with
#   the file "${comparison_file}". The executable is defined by
#
#     TEST_TARGET or
#     TEST_TARGET_DEBUG and TEST_TARGET_RELEASE
#
# - If the parameter file in the second test variant is named
#   "${test_name}.prm.in" it will be configured/preprocessed to a
#   "${test_name}.prm" file. This preprocessing is done with the CMake
#   macro CONFIGURE_FILE that replaces all strings @VARIABLE@ with the
#   contents of the corresponding CMake variable. This is useful in
#   particular to conveniently substitute @SOURCE_DIR@ with the full source
#   directory path of the test.
#
# For every deal.II build type (given by the variable DEAL_II_BUILD_TYPES)
# that is a (case insensitive) substring of CMAKE_BUILD_TYPE a test is
# defined.
#
# This macro gets the following options from the comparison file name (have
# a look at the testsuite documentation for details):
#  - usage of mpirun and number of simultaneous processes
#  - valid build configurations
#  - expected test stage
#
# This macro expects the CMAKE_BUILD_TYPE to be either "Debug", "Release",
# or "DebugRelease". For the first two build types, this macro will set up
# one test (linking against the corresponding deal.II library variant). In
# case of the build type "DebugRelease" two tests against both library
# variants will be set up. This macro throws a FATAL_ERROR if the deal.II
# installation does not support the requested library variant(s).
#
# The following variables must be set:
#
#   NUMDIFF_EXECUTABLE
#     - Complete path to the numdiff binary.
#
#   TEST_TIME_LIMIT
#     - specifying the maximal wall clock time in seconds a test is allowed
#       to run
#
# Usage:
#     DEAL_II_ADD_TEST(category test_name comparison_file)
#

MACRO(DEAL_II_ADD_TEST _category _test_name _comparison_file)

  IF(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    MESSAGE(FATAL_ERROR
      "\nDEAL_II_ADD_TEST can only be called in external (test sub-) projects after "
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

  #
  # Determine for which build types a test should be defined.
  #
  # Every deal.II build type (given by the list DEAL_II_BUILD_TYPES) that
  # is a (case insensitive) substring of CMAKE_BUILD_TYPE:
  #
  IF("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    SET(_build_types DEBUG)
  ELSEIF("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    SET(_build_types RELEASE)
  ELSEIF("${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    SET(_build_types DEBUG RELEASE)
  ELSE()
    MESSAGE(FATAL_ERROR
      "\nDEAL_II_ADD_TEST requires CMAKE_BUILD_TYPE to be set to "
      "\"Debug\", \"Release\", or \"DebugRelease\"\n\n"
      )
  ENDIF()

  FOREACH(_build ${_build_types})
    LIST(FIND DEAL_II_BUILD_TYPES ${_build} _match)
    IF("${_match}" STREQUAL "-1")
      MESSAGE(FATAL_ERROR
        "\nDEAL_II_ADD_TEST cannot set up a test with CMAKE_BUILD_TYPE "
        "\"${CMAKE_BUILD_TYPE}\". deal.II was build with CMAKE_BUILD_TYPE "
        "\"${DEAL_II_BUILD_TYPE}\"\n\n"
        )
    ENDIF()
  ENDFOREACH()

  FOREACH(_build ${_build_types})

    #
    # Obey "debug" and "release" keywords in the output file:
    #
    ITEM_MATCHES(_match "${_build}" ${_configuration})
    IF(_match OR "${_configuration}" STREQUAL "")

      STRING(TOLOWER ${_build} _build_lowercase)

      #
      # Select a suitable target:
      #
      IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.cc")

        SET(_source_file "${_test_name}.cc")

        SET(_target ${_test_name}.${_build_lowercase}) # target name
        SET(_run_args "$<TARGET_FILE:${_target}>") # the command to issue

      ELSEIF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.cu")

        SET(_source_file "${_test_name}.cu")

        SET(_target ${_test_name}.${_build_lowercase}) # target name
        SET(_run_args "$<TARGET_FILE:${_target}>") # the command to issue

      ELSEIF( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.prm" OR
              EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.prm.in" )

        IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.prm.in")
          SET(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
          CONFIGURE_FILE(
            "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.prm.in"
            "${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.prm"
            @ONLY
            )
          SET(_prm_file "${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.prm")
        ELSE()
          SET(_prm_file "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.prm")
        ENDIF()

        IF(NOT "${TEST_TARGET_${_build}}" STREQUAL "")
          SET(_target ${TEST_TARGET_${_build}})
        ELSEIF(NOT "${TEST_TARGET}" STREQUAL "")
          SET(_target ${TEST_TARGET})
        ELSE()
          MESSAGE(FATAL_ERROR
            "\nFor ${_comparison_file}: \"${_test_name}.prm(.in)\" provided, "
            "but neither \"\${TEST_TARGET}\", nor \"\${TEST_TARGET_${_build}}"
            "\" is defined.\n\n"
            )
        ENDIF()
        SET(_run_args
          "$<TARGET_FILE:${_target}>"
          "${_prm_file}"
          )

      ELSE()
        MESSAGE(FATAL_ERROR
          "\nFor ${_comparison_file}: Neither \"${_test_name}.cc\", "
          "nor \"${_test_name}.prm\" could be found!\n\n"
          )
      ENDIF()

      #
      # Set up a bunch of variables describing this particular test:
      #

      # If _n_cpu is equal to "0", a normal, sequential test will be run,
      # otherwise run the test with mpirun:
      IF("${_n_cpu}" STREQUAL "0")

        SET(_diff_target ${_test_name}.${_build_lowercase}.diff) # diff target name
        SET(_test_full ${_category}/${_test_name}.${_build_lowercase}) # full test name
        SET(_test_directory ${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.${_build_lowercase}) # directory to run the test in

      ELSE()

        SET(_diff_target ${_test_name}.mpirun${_n_cpu}.${_build_lowercase}.diff) # diff target name
        SET(_test_full ${_category}/${_test_name}.mpirun=${_n_cpu}.${_build_lowercase}) # full test name
        SET(_test_directory ${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.${_build_lowercase}/mpirun=${_n_cpu}) # directory to run the test in
        SET(_run_args
          "${DEAL_II_MPIEXEC}"
          ${DEAL_II_MPIEXEC_NUMPROC_FLAG} ${_n_cpu}
          ${DEAL_II_MPIEXEC_PREFLAGS}
          ${_run_args}
          "${DEAL_II_MPIEXEC_POSTFLAGS}"
          )
      ENDIF()

      FILE(MAKE_DIRECTORY ${_test_directory})

      #
      # Add an executable (for the first type of tests) and set up compile
      # definitions and the full link interface. Only add the target once.
      #

      IF(NOT TARGET ${_target})
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


        ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL
          ${_generated_files}
          ${_source_file}
          ${CMAKE_CURRENT_BINARY_DIR}/${_target}/interrupt_guard.cc
          )

        DEAL_II_SETUP_TARGET(${_target} ${_build})
        TARGET_LINK_LIBRARIES(${_target}
          ${TEST_LIBRARIES} ${TEST_LIBRARIES_${_build}}
          )

        SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
          COMPILE_DEFINITIONS SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
          )
        SET_PROPERTY(TARGET ${_target} PROPERTY
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_target}"
          )

      ENDIF()

      #
      # Add a top level target to run and compare the test:
      #

      ADD_CUSTOM_COMMAND(OUTPUT ${_test_directory}/output
        COMMAND sh ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.sh
          run "${_test_full}" ${_run_args}
        COMMAND ${PERL_EXECUTABLE}
          -pi ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/normalize.pl
          ${_test_directory}/output
        WORKING_DIRECTORY
          ${_test_directory}
        DEPENDS
          ${_target}
          ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/normalize.pl
        VERBATIM
        )

      FILE(GLOB _comparison_files ${_comparison_file} ${_comparison_file}.*)

      ADD_CUSTOM_COMMAND(OUTPUT ${_test_directory}/diff
        COMMAND sh ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.sh
          diff "${_test_full}" "${NUMDIFF_EXECUTABLE}"
          "${_comparison_file}" ${_run_args}
        WORKING_DIRECTORY
          ${_test_directory}
        DEPENDS
          ${_test_directory}/output
          ${_comparison_files}
        VERBATIM
        )

      ADD_CUSTOM_TARGET(${_diff_target}
        COMMAND echo "${_test_full}: BUILD successful."
        COMMAND echo "${_test_full}: RUN successful."
        COMMAND echo "${_test_full}: DIFF successful."
        COMMAND echo "${_test_full}: PASSED."
        DEPENDS ${_test_directory}/diff
        )

      #
      # And finally define the test:
      #

      ADD_TEST(NAME ${_test_full}
        COMMAND ${CMAKE_COMMAND}
          -DTRGT=${_diff_target}
          -DTEST=${_test_full}
          -DEXPECT=${_expect}
          -DBINARY_DIR=${CMAKE_BINARY_DIR}
          -DGUARD_FILE=${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.${_build_lowercase}/interrupt_guard.cc
          -P ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.cmake
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
