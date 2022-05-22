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
#     - Specifies the maximal wall clock time in seconds a test is allowed
#       to run.
#
#   TEST_MPI_RANK_LIMIT
#     - Specifies the maximal number of MPI ranks that can be used. If a
#       test variant configures a larger number of MPI ranks (via
#       .mpirun=N. in the output file) than this limit the test will be
#       dropped. The special value 0 enforces no limit. Defaults to 0.
#
#   TEST_THREAD_LIMIT
#     - Specifies the maximal number of worker threads that can should be
#       used by the threading backend. If a test variant configures a
#       larger number of threads (via .threads=N. in the output file) than
#       this limit the test will be dropped. Note that individual tests
#       might exceed this limit by calling
#       MultithreadInfo::set_thread_limit(), or by manually creating
#       additional threads. The special value 0 enforces no limit. Defaults
#       to 0.
#
#   ENABLE_PERFORMANCE_TESTS
#     - If defined and set to true the execution of performance tests will
#       be enabled.
#
#   TESTING_ENVIRONMENT
#     - Specifies the performance test testing environment. Valid options
#       are:
#         * "light":  mobile laptop, >=2 physical cores, >=8GB RAM
#         * "medium": workstation, >=8 physical cores, >=32GB RAM
#         * "heavy":  compute node, >=32 physical cores, >=128GB RAM
#
# Usage:
#     DEAL_II_ADD_TEST(category test_name comparison_file)
#

FUNCTION(DEAL_II_ADD_TEST _category _test_name _comparison_file)

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
  STRING(REGEX MATCH "mpirun=([0-9]+|max)" _n_cpu ${_file})
  IF("${_n_cpu}" STREQUAL "")
    SET(_n_cpu 0) # 0 indicates that no mpirun should be used
  ELSE()
    STRING(REGEX REPLACE "^mpirun=([0-9]+|max)$" "\\1" _n_cpu ${_n_cpu})
  ENDIF()

  #
  # If we encounter the special string "mpirun=max" set the number of MPI
  # ranks used for the test to the maximum number of allowed ranks. If no
  # limit has been specified, i.e., TEST_MPI_RANK_LIMIT is 0, skip defining
  # the test.
  #
  IF("${_n_cpu}" STREQUAL "max")
    IF(TEST_MPI_RANK_LIMIT EQUAL 0)
      RETURN()
    ENDIF()
    SET(_n_cpu "${TEST_MPI_RANK_LIMIT}")
  ENDIF()

  #
  # If the number of MPI ranks specified for the test via .mpirun=N.
  # exceeds the limit ${TEST_MPI_RANK_LIMIT}, skip defining the test
  #
  IF(TEST_MPI_RANK_LIMIT GREATER 0 AND _n_cpu GREATER TEST_MPI_RANK_LIMIT)
    RETURN()
  ENDIF()

  #
  # Determine whether the test declaration specifies a thread pool size via
  # threads=N:
  #
  STRING(REGEX MATCH "threads=([0-9]+|max)" _n_threads ${_file})
  IF("${_n_threads}" STREQUAL "")
    SET(_n_threads 0) # 0 indicates that the default thread pool size
                      # should be used (currently set to 3 in tests.h)
  ELSE()
    STRING(REGEX REPLACE "^threads=([0-9]+|max)$" "\\1" _n_threads ${_n_threads})
  ENDIF()

  #
  # If we encounter the special string "threads=max" set the number of
  # threads of the threading pool to the maximum number of allowed threads.
  # If no limit has been specified, i.e., TEST_THREAD_LIMIT is 0, skip
  # defining the test.
  #
  IF("${_n_threads}" STREQUAL "max")
    IF(TEST_THREAD_LIMIT EQUAL 0)
      RETURN()
    ENDIF()
    SET(_n_threads "${TEST_THREAD_LIMIT}")
  ENDIF()

  #
  # If the number of threads specified for the test via .threads=N. exceeds
  # the limit ${TEST_THREAD_LIMIT}, skip defining the test
  #
  IF(TEST_THREAD_LIMIT GREATER 0 AND _n_threads GREATER TEST_THREAD_LIMIT)
    RETURN()
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
  # Determine whether the .run_only keyword is present:
  #
  SET(_run_only FALSE)
  IF(_file MATCHES "\\.run_only$")
    SET(_run_only TRUE)
  ENDIF()

  #
  # Determine whether the .exclusive. keyword is present:
  #
  SET(_exclusive FALSE)
  IF(_file MATCHES "\\.exclusive\\.")
    SET(_exclusive TRUE)
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

  IF(NOT TESTING_ENVIRONMENT MATCHES "^(light|medium|heavy)$")
    MESSAGE(FATAL_ERROR
      "The TESTING_ENVIRONMENT variable must be set to either \"light\","
      " \"medium\", or \"heavy\"."
      )
  ENDIF()

  #
  # Determine source or parameter file:
  #

  FILE(GLOB _source_file
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.c[cu]"
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.prm"
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.prm.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.json"
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.json.in"
    )

  LIST(LENGTH _source_file _number)
  IF(NOT _number EQUAL 1)
    IF(_number EQUAL 0)
      MESSAGE(FATAL_ERROR "\n${_comparison_file}:\n"
        "A comparison file (ending in .output or .run-only) has been "
        "picked up but no suitable source file or parameter file was "
        "found. Please provide exactly one of the following.\n"
        "A source file \"${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.c[cu]\",\n"
        "or a parameter file \"${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.(prm|json)[.in]\".\n"
        )
    ELSE()
      STRING(REPLACE ";" "\n" _source_file "${_source_file}")
      MESSAGE(FATAL_ERROR "\n${_comparison_file}:\n"
        "A comparison file (ending in .output or .run-only) has been "
        "picked up with multiple suitable source files or parameter files:\n"
        "${_source_file}\n"
        "There must be exactly one source or parameter file.\n"
        )
    ENDIF()
  ENDIF()

  #
  # Run CONFIGURE_FILE on every parameter file ending in "in"
  #

  IF("${_source_file}" MATCHES ".in$")
    SET(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    STRING(REGEX MATCH "(json|prm).in$" _suffix "${_source_file}")
    STRING(REPLACE ".in" "" _suffix "${_suffix}")
    CONFIGURE_FILE(
      "${_source_file}"
      "${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.${_suffix}"
      @ONLY
      )
    SET(_source_file "${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.${_suffix}")
  ENDIF()

  FOREACH(_build ${_build_types})
    #
    # Obey "debug" and "release" keywords in the output file:
    #
    ITEM_MATCHES(_match "${_build}" ${_configuration})
    IF(_match OR "${_configuration}" STREQUAL "")
      STRING(TOLOWER ${_build} _build_lowercase)

      SET(_target ${_category}.${_test_name}.${_build_lowercase}) # target name
      SET(_target_short ${_test_name}.${_build_lowercase}) # short target name
      SET(_run_args "$<TARGET_FILE:${_target}>") # the command to issue

      #
      # If the variable ${category_test_RUNARGS_PREFIX) is nonempty prepend
      # it to the command line of the test:
      #
      if(NOT "${${_category}_${_test_name}_RUNARGS_PREFIX}" STREQUAL "")
        SET(_run_args
          ${${_category}_${_test_name}_RUNARGS_PREFIX} "${_run_args}")
      ENDIF()

      #
      # Override target and run command for parameter file variants:
      #
      IF("${_source_file}" MATCHES "(prm|json)$")
        IF(NOT "${TEST_TARGET_${_build}}" STREQUAL "")
          SET(_target ${TEST_TARGET_${_build}})
        ELSEIF(NOT "${TEST_TARGET}" STREQUAL "")
          SET(_target ${TEST_TARGET})
        ELSE()
          MESSAGE(FATAL_ERROR "\n${_comparison_file}:\n"
            "A parameter file \"${_test_name}.(prm|json)(|.in)\" has been "
            "found, but neither \"\${TEST_TARGET}\", nor "
            "\"\${TEST_TARGET_${_build}}\" have been defined.\n"
            )
        ENDIF()
        SET(_target_short ${_target})
        SET(_run_args
          "$<TARGET_FILE:${_target}>"
          "${_source_file}"
          )
      ENDIF()

      #
      # If _n_cpu or _n_threads are larger than zero we have to accomodate
      # the fact that multiple output files specifying a different mpirun or
      # threads count are present. In order to accomodate this we create a
      # runtime subdirectory "mpirun_M-threads_N" to the test.
      #
      # Note that we could do this unconditionally for every test but for
      # aesthetic reasons chose to not create the directory for
      # "mpirun_0-threads_0".
      #

      SET(_test_target    ${_category}.${_test_name}) # diff target name
      SET(_test_full      ${_category}/${_test_name}) # full test name
      SET(_test_directory ${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.${_build_lowercase}) # directory to run the test in

      IF(NOT "${_n_cpu}" STREQUAL "0")
        STRING(APPEND _test_target   ".mpirun${_n_cpu}")
        STRING(APPEND _test_full     ".mpirun=${_n_cpu}")
        STRING(APPEND _test_directory "/mpirun=${_n_cpu}")
      ENDIF()

      IF(NOT "${_n_threads}" STREQUAL "0")
        STRING(APPEND _test_target   ".threads${_n_threads}")
        STRING(APPEND _test_full     ".threads=${_n_threads}")
        STRING(APPEND _test_directory "/threads=${_n_threads}")
      ENDIF()

      STRING(APPEND _test_target ".${_build_lowercase}.test")
      STRING(APPEND _test_full   ".${_build_lowercase}")

      #
      # Test variants with ".mpirun=[...]." have to be executed via mpirun
      # (or whatever ${DEAL_II_MPIEXEC} is set to).
      #
      IF(NOT "${_n_cpu}" STREQUAL "0")
        SET(_run_args
          "${DEAL_II_MPIEXEC}"
          ${DEAL_II_MPIEXEC_NUMPROC_FLAG} ${_n_cpu}
          ${DEAL_II_MPIEXEC_PREFLAGS}
          ${_run_args}
          ${DEAL_II_MPIEXEC_POSTFLAGS}
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
          OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_target_short}/interrupt_guard.cc
          COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/${_target_short}/interrupt_guard.cc
          )

        ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL
          ${_generated_files}
          ${_source_file}
          ${CMAKE_CURRENT_BINARY_DIR}/${_target_short}/interrupt_guard.cc
          )

        SET_TARGET_PROPERTIES(${_target} PROPERTIES OUTPUT_NAME ${_target_short})

        DEAL_II_SETUP_TARGET(${_target} ${_build})
        TARGET_LINK_LIBRARIES(${_target}
          ${TEST_LIBRARIES} ${TEST_LIBRARIES_${_build}}
          )

        SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
          COMPILE_DEFINITIONS
            SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
            TESTING_ENVIRONMENT=${TESTING_ENVIRONMENT}
          )

        IF(ENABLE_PERFORMANCE_TESTS)
          SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
            COMPILE_DEFINITIONS ENABLE_PERFORMANCE_TESTS
            )
        ENDIF()

        SET_PROPERTY(TARGET ${_target} PROPERTY
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_target_short}"
          )

      ENDIF()

      #
      # Add a top level target to run and compare the test:
      #

      ADD_CUSTOM_COMMAND(OUTPUT ${_test_directory}/output
        COMMAND TEST_N_THREADS=${_n_threads}
          sh ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.sh
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

      IF(_run_only)
        ADD_CUSTOM_TARGET(${_test_target}
          COMMAND echo "${_test_full}: BUILD successful."
          COMMAND echo "${_test_full}: RUN successful."
          COMMAND echo "${_test_full}: DIFF skipped."
          COMMAND echo "${_test_full}: PASSED."
          DEPENDS ${_test_directory}/output
          )

      ELSE()

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

        ADD_CUSTOM_TARGET(${_test_target}
          COMMAND echo "${_test_full}: BUILD successful."
          COMMAND echo "${_test_full}: RUN successful."
          COMMAND echo "${_test_full}: DIFF successful."
          COMMAND echo "${_test_full}: PASSED."
          DEPENDS ${_test_directory}/diff
          )

      ENDIF()

      #
      # And finally define the test:
      #

      ADD_TEST(NAME ${_test_full}
        COMMAND ${CMAKE_COMMAND}
          -DTRGT=${_test_target}
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

      IF(_exclusive)
        #
        # Ensure that the test is not executed concurrently with any other
        # tests.
        #
        SET_TESTS_PROPERTIES(${_test_full} PROPERTIES RUN_SERIAL TRUE)

      ELSEIF(NOT ENABLE_PERFORMANCE_TESTS)
        #
        # Limit concurrency of mpi tests. We can only set concurrency for
        # the entire test, which includes the compiling and linking stages
        # that are purely sequential. There is no good way to model this
        # without unnecessarily restricting concurrency. Consequently, we
        # just choose to model an "average" concurrency as one half of the
        # number of MPI jobs.
        #
        IF(_n_cpu GREATER 2)
          MATH(EXPR _slots "${_n_cpu} / 2")
          SET_TESTS_PROPERTIES(${_test_full} PROPERTIES PROCESSORS ${_slots})
        ENDIF()

      ELSE()
        #
        # In case ENABLE_PERFORMANCE_TESTS is set we limit the concurrency
        # of performance tests to the number of specified mpi ranks times
        # the number of specified threads.
        #
        SET(_slots 1)
        IF(_n_cpu GREATER 0)
          MATH(EXPR _slots "${_slots} * ${_n_cpu}")
        ENDIF()
        IF(_n_threads GREATER 0)
          MATH(EXPR _slots "${_slots} * ${_n_threads}")
        ENDIF()
        SET_TESTS_PROPERTIES(${_test_full} PROPERTIES PROCESSORS ${_slots})
      ENDIF()

      IF(NOT "${_n_cpu}${_n_threads}" STREQUAL "00")
        #
        # Running multiple variants in parallel triggers a race condition
        # where the same (not yet existent) executable is built
        # concurrently leading to undefined outcomes.
        #
        # Luckily CMake has a mechanism to force a test to be run after
        # another has finished (and both are scheduled):
        #
        IF(DEFINED TEST_DEPENDENCIES_${_target})
          SET_TESTS_PROPERTIES(${_test_full} PROPERTIES
            DEPENDS ${TEST_DEPENDENCIES_${_target}}
            )
        ENDIF()
        SET(TEST_DEPENDENCIES_${_target} ${_test_full} PARENT_SCOPE)
      ENDIF()

    ENDIF()
  ENDFOREACH()
ENDFUNCTION()
