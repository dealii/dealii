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
#   BASH
#     - Complete path to the bash shell.
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
#     - If set to true the execution of performance tests will be enabled.
#
#   TESTING_ENVIRONMENT
#     - Specifies the performance test testing environment. Valid options
#       are:
#         * "light":  mobile laptop, >=2 physical cores, >=8GB RAM
#         * "medium": workstation, >=8 physical cores, >=32GB RAM
#         * "heavy":  compute node, >=32 physical cores, >=128GB RAM
#
# Usage:
#     deal_ii_add_test(category test_name comparison_file)
#

#
# A small helper macro that is used below:
#

macro(item_matches _var _regex)
  set(${_var})
  foreach (_item ${ARGN})
    if("${_item}" MATCHES ${_regex})
      set(${_var} TRUE)
      break()
    endif()
  endforeach()
endmacro()


function(deal_ii_add_test _category _test_name _comparison_file)
  if(NOT TARGET compile_test_executables)
    add_custom_target(compile_test_executables)
  endif()

  if(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    message(FATAL_ERROR
      "\nDEAL_II_ADD_TEST can only be called in external (test sub-) projects after "
      "the inclusion of deal.IIConfig.cmake. It is not intended for "
      "internal use.\n\n"
      )
  endif()

  get_filename_component(_file ${_comparison_file} NAME)

  #
  # Determine valid build configurations for this test:
  #
  set(_configuration)
  if(_file MATCHES "\\.debug\\.")
    set(_configuration DEBUG)
  elseif(_file MATCHES "\\.release\\.")
    set(_configuration RELEASE)
  endif()

  #
  # Determine whether the test should be run with mpirun:
  #
  string(REGEX MATCH "mpirun=([0-9]+|max)" _n_cpu ${_file})
  if("${_n_cpu}" STREQUAL "")
    set(_n_cpu 0) # 0 indicates that no mpirun should be used
  else()
    string(REGEX REPLACE "^mpirun=([0-9]+|max)$" "\\1" _n_cpu ${_n_cpu})
  endif()

  #
  # If we encounter the special string "mpirun=max" set the number of MPI
  # ranks used for the test to the maximum number of allowed ranks. If no
  # limit has been specified, i.e., TEST_MPI_RANK_LIMIT is 0, skip defining
  # the test.
  #
  # This mechanism is specifically used in the performance test suite.
  #
  if("${_n_cpu}" STREQUAL "max")
    if(TEST_MPI_RANK_LIMIT EQUAL 0)
      return()
    endif()
    set(_n_cpu "${TEST_MPI_RANK_LIMIT}")
  endif()

  #
  # If the number of MPI ranks specified for the test via .mpirun=N.
  # exceeds the limit ${TEST_MPI_RANK_LIMIT}, skip defining the test
  #
  if(TEST_MPI_RANK_LIMIT GREATER 0 AND _n_cpu GREATER TEST_MPI_RANK_LIMIT)
    return()
  endif()

  #
  # Determine whether the test declaration specifies a thread pool size via
  # threads=N:
  #
  string(REGEX MATCH "threads=([0-9]+|max)" _n_threads ${_file})
  if("${_n_threads}" STREQUAL "")
    set(_n_threads 0) # 0 indicates that the default thread pool size
                      # should be used (currently set to 3 in tests.h)
  else()
    string(REGEX REPLACE "^threads=([0-9]+|max)$" "\\1" _n_threads ${_n_threads})
  endif()

  #
  # If we encounter the special string "threads=max" set the number of
  # threads of the threading pool to the maximum number of allowed threads.
  # If no limit has been specified, i.e., TEST_THREAD_LIMIT is 0, skip
  # defining the test.
  #
  if("${_n_threads}" STREQUAL "max")
    if(TEST_THREAD_LIMIT EQUAL 0)
      return()
    endif()
    set(_n_threads "${TEST_THREAD_LIMIT}")
  endif()

  #
  # If the number of threads specified for the test via .threads=N. exceeds
  # the limit ${TEST_THREAD_LIMIT}, skip defining the test
  #
  if(TEST_THREAD_LIMIT GREATER 0 AND _n_threads GREATER TEST_THREAD_LIMIT)
    return()
  endif()

  #
  # Determine whether the .run_only keyword is present.
  #
  # In case no numdiff executable was found we fall back to simply running
  # the tests as well (but not comparing them).
  #
  set(_run_only FALSE)
  if(_file MATCHES "\\.run_only$" OR "${NUMDIFF_EXECUTABLE}" STREQUAL "")
    set(_run_only TRUE)
  endif()

  #
  # Determine the expected build stage of this test:
  #
  string(REGEX MATCH "expect=([a-z]*)" _expect ${_file})
  if("${_expect}" STREQUAL "")
    set(_expect "PASSED")
  else()
    string(REGEX REPLACE "^expect=([a-z]*)$" "\\1" _expect ${_expect})
    string(TOUPPER ${_expect} _expect)
  endif()

  #
  # If _run_only is set then we don't compare test results. Therefore,
  # we won't fail in the "DIFF" stage of a test.
  #
  if(_run_only AND "${_expect}" STREQUAL "DIFF")
    set(_expect "PASSED")
  endif()

  #
  # Determine whether the .exclusive. keyword is present:
  #
  set(_exclusive FALSE)
  if(_file MATCHES "\\.exclusive\\.")
    set(_exclusive TRUE)
  endif()

  #
  # Determine for which build types a test should be defined.
  #
  # Every deal.II build type (given by the list DEAL_II_BUILD_TYPES) that
  # is a (case insensitive) substring of CMAKE_BUILD_TYPE:
  #
  if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(_build_types DEBUG)
  elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    set(_build_types RELEASE)
  elseif("${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    set(_build_types DEBUG RELEASE)
  else()
    message(FATAL_ERROR
      "\nDEAL_II_ADD_TEST requires CMAKE_BUILD_TYPE to be set to "
      "\"Debug\", \"Release\", or \"DebugRelease\"\n\n"
      )
  endif()

  foreach(_build ${_build_types})
    list(FIND DEAL_II_BUILD_TYPES ${_build} _match)
    if("${_match}" STREQUAL "-1")
      message(FATAL_ERROR
        "\nDEAL_II_ADD_TEST cannot set up a test with CMAKE_BUILD_TYPE "
        "\"${CMAKE_BUILD_TYPE}\". deal.II was build with CMAKE_BUILD_TYPE "
        "\"${DEAL_II_BUILD_TYPE}\"\n\n"
        )
    endif()
  endforeach()

  if(NOT TESTING_ENVIRONMENT MATCHES "^(light|medium|heavy)$")
    message(FATAL_ERROR
      "The TESTING_ENVIRONMENT variable must be set to either \"light\","
      " \"medium\", or \"heavy\"."
      )
  endif()

  #
  # Determine source or parameter file:
  #

  file(GLOB _source_file
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.c[cu]"
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.prm"
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.prm.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.json"
    "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.json.in"
    )

  list(LENGTH _source_file _number)
  if(NOT _number EQUAL 1)
    if(_number EQUAL 0)
      #
      # None of the candidates above exist. Check whether a
      # ${_test_name}.cc file gets generated during the build process:
      set(_source_file "${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.cc")
      get_property(_generated SOURCE "${_source_file}" PROPERTY GENERATED)
      if(NOT ${_generated})
        message(FATAL_ERROR "\n${_comparison_file}:\n"
          "A comparison file (ending in .output or .run-only) has been "
          "picked up but no suitable source file, generated source file, or "
          "parameter file was found. Please provide exactly one of the following.\n"
          "A source file \"${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.c[cu]\",\n"
          "or a parameter file \"${CMAKE_CURRENT_SOURCE_DIR}/${_test_name}.(prm|json)[.in]\".\n"
          )
      endif()

    else()
      string(REPLACE ";" "\n" _source_file "${_source_file}")
      message(FATAL_ERROR "\n${_comparison_file}:\n"
        "A comparison file (ending in .output or .run-only) has been "
        "picked up with multiple suitable source files or parameter files:\n"
        "${_source_file}\n"
        "There must be exactly one source or parameter file.\n"
        )
    endif()
  endif()

  #
  # Run CONFIGURE_FILE on every parameter file ending in "in"
  #

  if("${_source_file}" MATCHES ".in$")
    set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    string(REGEX MATCH "(json|prm).in$" _suffix "${_source_file}")
    string(REPLACE ".in" "" _suffix "${_suffix}")
    configure_file(
      "${_source_file}"
      "${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.${_suffix}"
      @ONLY
      )
    set(_source_file "${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.${_suffix}")
  endif()

  foreach(_build ${_build_types})
    #
    # Increment the _number_of_tests counter by one in the parent scope.
    # Note that the variable can be undefined initially. In this case we
    # simply store "1" in the first iteration.
    #
    math(EXPR _number_of_tests "${_number_of_tests} + 1")
    set(_number_of_tests "${_number_of_tests}" PARENT_SCOPE)

    #
    # Obey "debug" and "release" keywords in the output file:
    #
    item_matches(_match "${_build}" ${_configuration})
    if(_match OR "${_configuration}" STREQUAL "")
      string(TOLOWER ${_build} _build_lowercase)

      set(_target ${_category}.${_test_name}.${_build_lowercase}) # target name
      set(_target_short ${_test_name}.${_build_lowercase}) # short target name
      set(_run_args "$<TARGET_FILE:${_target}>") # the command to issue

      #
      # If the variable ${category_test_RUNARGS_PREFIX) is nonempty prepend
      # it to the command line of the test:
      #
      if(NOT "${${_category}_${_test_name}_RUNARGS_PREFIX}" STREQUAL "")
        set(_run_args
          ${${_category}_${_test_name}_RUNARGS_PREFIX} "${_run_args}")
      endif()

      #
      # Override target and run command for parameter file variants:
      #
      if("${_source_file}" MATCHES "(prm|json)$")
        if(TARGET "${TEST_TARGET_${_build}}")
          set(_target ${TEST_TARGET_${_build}})
        elseif(TARGET "${TEST_TARGET}")
          set(_target ${TEST_TARGET})
        else()
          message(FATAL_ERROR
            "The parameter file \"${_source_file}\" has been found, but neither "
            "\"\${TEST_TARGET}\", nor \"\${TEST_TARGET_${_build}}\" have been set to "
            "valid target name.\n"
            )
        endif()
        set(_target_short ${_target})
        set(_run_args
          "$<TARGET_FILE:${_target}>"
          "${_source_file}"
          )
      endif()

      #
      # If _n_cpu or _n_threads are larger than zero we have to accommodate
      # the fact that multiple output files specifying a different mpirun or
      # threads count are present. In order to accommodate this we create a
      # runtime subdirectory "mpirun_M-threads_N" to the test.
      #
      # Note that we could do this unconditionally for every test but for
      # aesthetic reasons chose to not create the directory for
      # "mpirun_0-threads_0".
      #

      set(_test_target    ${_category}.${_test_name}) # diff/run target name
      set(_test_full      ${_category}/${_test_name}) # full test name
      set(_test_directory ${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.${_build_lowercase}) # directory to run the test in

      if("${_n_cpu}" STREQUAL "0" AND "${_n_threads}" STREQUAL "0")
        string(APPEND _test_directory "/serial")
      endif()

      if(NOT "${_n_cpu}" STREQUAL "0")
        string(APPEND _test_target   ".mpirun${_n_cpu}")
        string(APPEND _test_full     ".mpirun=${_n_cpu}")
        string(APPEND _test_directory "/mpirun=${_n_cpu}")
      endif()

      if(NOT "${_n_threads}" STREQUAL "0")
        string(APPEND _test_target   ".threads${_n_threads}")
        string(APPEND _test_full     ".threads=${_n_threads}")
        string(APPEND _test_directory "/threads=${_n_threads}")
      endif()

      string(APPEND _test_target ".${_build_lowercase}.test")
      string(APPEND _test_full   ".${_build_lowercase}")

      #
      # Test variants with ".mpirun=[...]." have to be executed via mpirun
      # (or whatever ${DEAL_II_MPIEXEC} is set to).
      #
      if(NOT "${_n_cpu}" STREQUAL "0")
        set(_run_args
          "${DEAL_II_MPIEXEC}"
          ${DEAL_II_MPIEXEC_NUMPROC_FLAG} ${_n_cpu}
          ${DEAL_II_MPIEXEC_PREFLAGS}
          ${_run_args}
          ${DEAL_II_MPIEXEC_POSTFLAGS}
          )
      endif()

      file(MAKE_DIRECTORY ${_test_directory})

      #
      # Determine whether the test shares a common executable target. This
      # involves tests with .threads=N. and .mpirun=N. annotation, as well
      # as tests with parameter files (that might share a common executable
      # target).
      #
      # In this case we have to make sure that concurrently invoking the
      # test does not accidentally trigger a concurrent build of the
      # executable target. We ensure this by declaring an additional test
      # that only builds the shared target / ensures the shared target is
      # present. All run tests then requires this test target as a "setup
      # fixture", see
      # https://cmake.org/cmake/help/latest/prop_test/FIXTURES_REQUIRED.html#prop_test:FIXTURES_REQUIRED
      #
      set(_shared_target FALSE)
      if(NOT "${_n_cpu}${_n_threads}" STREQUAL "00" OR "${_source_file}" MATCHES "(prm|json)$")
        set(_shared_target TRUE)

        #
        # Build system-internal target name and final test name for the
        # "executable" test. We have to make sure that the target and test
        # names stay the same independent of test name and test category,
        # thus the rather funny name:
        #
        set(_test_executable_target "test_dependency.${_target}.executable")
        set(_test_executable_full   "test_dependency/${_target}.executable")
      endif()

      #
      # Add an executable (for the first type of tests) and set up compile
      # definitions and the full link interface. Only add the target once.
      #

      if(NOT TARGET ${_target})

        add_executable(${_target} EXCLUDE_FROM_ALL
          ${_generated_files}
          ${_source_file}
          )

        set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${_target_short})

        deal_ii_setup_target(${_target} ${_build})
        target_link_libraries(${_target}
          ${TEST_LIBRARIES} ${TEST_LIBRARIES_${_build}}
          )

        set_property(TARGET ${_target} APPEND PROPERTY
          COMPILE_DEFINITIONS
            SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
            TESTING_ENVIRONMENT=${TESTING_ENVIRONMENT}
          )

        if(ENABLE_PERFORMANCE_TESTS)
          set_property(TARGET ${_target} APPEND PROPERTY
            COMPILE_DEFINITIONS ENABLE_PERFORMANCE_TESTS
            )
        endif()

        set_property(TARGET ${_target} PROPERTY
          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_target_short}"
          )

        add_dependencies(compile_test_executables ${_target})
      endif()

      #
      # Add a top level target to compile the test:
      #

      if(_shared_target AND NOT TARGET ${_test_executable_target})
        #
        # Increment the _number_of_test_dependencies counter by one in the
        # parent scope. Test dependencies are all tests where we have split
        # out compiling and linking of an executable into a separate
        # "test_dependency/" test. Note that the variable can be undefined
        # initially. In this case we simply store "1" in the first
        # iteration.
        #
        math(EXPR _number_of_test_dependencies "${_number_of_test_dependencies} + 1")
        set(_number_of_test_dependencies "${_number_of_test_dependencies}" PARENT_SCOPE)

        add_custom_target(${_test_executable_target}
          COMMAND echo "${_test_executable_full}: BUILD successful."
          COMMAND echo "${_test_executable_full}: RUN skipped."
          COMMAND echo "${_test_executable_full}: DIFF skipped."
          COMMAND echo "${_test_executable_full}: PASSED."
          DEPENDS ${_target}
          )
        add_test(NAME ${_test_executable_full}
          COMMAND ${CMAKE_COMMAND}
            -DTRGT=${_test_executable_target}
            -DTEST=${_test_executable_full}
            -DEXPECT=PASSED
            -DBINARY_DIR=${CMAKE_BINARY_DIR}
            -P ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.cmake
          WORKING_DIRECTORY ${_test_directory}
          )
        set_tests_properties(${_test_executable_full} PROPERTIES
          LABEL "test_dependency"
          TIMEOUT ${TEST_TIME_LIMIT}
          FIXTURES_SETUP ${_test_executable_full}
          )
      endif()

      #
      # Add a top level target to run and compare the test:
      #

      add_custom_command(OUTPUT ${_test_directory}/output
        COMMAND ${BASH} ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.sh
          run "${_test_full}" TEST_N_THREADS=${_n_threads} ${_run_args}
        COMMAND ${PERL_EXECUTABLE}
          -pi ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/normalize.pl
          ${_test_directory}/output
        WORKING_DIRECTORY
          ${_test_directory}
        DEPENDS
          ${_target}
          ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/normalize.pl
        COMMENT
          "Normalizing test output file ${_test_directory}/output"
        VERBATIM
        )

      if(_run_only)
        #
        # Only compile and run the test executable. Do not run a diff
        # stage. We use this feature for performance tests (where comparing
        # output does not make sense), or for tests that signal success or
        # failure with a return code (such as our quick tests).
        #

        add_custom_target(${_test_target}
          COMMAND echo "${_test_full}: BUILD successful."
          COMMAND echo "${_test_full}: RUN successful."
          COMMAND echo "${_test_full}: DIFF skipped."
          COMMAND echo "${_test_full}: PASSED."
          DEPENDS ${_test_directory}/output
          )

      else()
        #
        # Add a diff rule and set up a test target that depends on a
        # successful compilation, run and diff.
        #

        file(GLOB _comparison_files ${_comparison_file} ${_comparison_file}.*)

        add_custom_command(OUTPUT ${_test_directory}/diff
          COMMAND ${BASH} ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.sh
            diff "${_test_full}" "${NUMDIFF_EXECUTABLE}"
            "${_comparison_file}" ${_run_args}
          WORKING_DIRECTORY
            ${_test_directory}
          DEPENDS
            ${_test_directory}/output
            ${_comparison_files}
          VERBATIM
          )

        add_custom_target(${_test_target}
          COMMAND echo "${_test_full}: BUILD successful."
          COMMAND echo "${_test_full}: RUN successful."
          COMMAND echo "${_test_full}: DIFF successful."
          COMMAND echo "${_test_full}: PASSED."
          DEPENDS ${_test_directory}/diff
          )

      endif()

      #
      # And finally define the test:
      #

      add_test(NAME ${_test_full}
        COMMAND ${CMAKE_COMMAND}
          -DTRGT=${_test_target}
          -DTEST=${_test_full}
          -DEXPECT=${_expect}
          -DBINARY_DIR=${CMAKE_BINARY_DIR}
          -P ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.cmake
        WORKING_DIRECTORY ${_test_directory}
        )
      set_tests_properties(${_test_full} PROPERTIES
        LABEL "${_category}"
        TIMEOUT ${TEST_TIME_LIMIT}
        )

      if(_shared_target)
        set_tests_properties(${_test_full} PROPERTIES
          FIXTURES_REQUIRED ${_test_executable_full}
          )
      endif()

      if(_exclusive)
        #
        # Ensure that the test is not executed concurrently with any other
        # tests.
        #
        set_tests_properties(${_test_full} PROPERTIES
          RUN_SERIAL TRUE
          ENVIRONMENT "TEST_IS_EXCLUSIVE=true"
          )

      else()
        #
        # Limit concurrency of tests that run on multiple mpi ranks, or
        # that explicitly spawn multiple worker threads.
        #
        set(_slots 1)
        if(_n_cpu GREATER 0)
          math(EXPR _slots "${_slots} * ${_n_cpu}")
        endif()
        if(_n_threads GREATER 0)
          math(EXPR _slots "${_slots} * ${_n_threads}")
        endif()
        set_tests_properties(${_test_full} PROPERTIES PROCESSORS ${_slots})
      endif()
    endif()
  endforeach()
endfunction()
