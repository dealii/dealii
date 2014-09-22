## ---------------------------------------------------------------------
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
INCLUDE(${DEAL_II_SOURCE_DIR}/tests/setup_testsubproject.cmake)
PROJECT(testsuite CXX)
INCLUDE(${DEAL_II_TARGET_CONFIG})

#
# build_tests tests are special:
#
# Try to configure and build every example step in debug and release
# configuration. Error condition is that if a test is successfully
# configured it must compile successfully...
#
# ... except for the tests listed below: They have to configure, build and
# run successfully in the given build configuration:
#

SET(_category build_tests)

# Run a minimalistic set of steps in debug configuration:
SET(_debug_steps
  step-1  step-2  step-3  step-4  step-5
  step-6  step-47
  )

# Run all configuration independent steps in release configuration:
SET(_release_steps
  step-1  step-2  step-3  step-4  step-5
  step-6  step-7  step-8  step-9  step-10
  step-11 step-12 step-13 step-14 step-16
  step-20 step-23 step-25 step-26 step-27
  step-30 step-38 step-39 step-44 step-45
  step-47 step-48 step-49
  )

FILE(GLOB _steps ${DEAL_II_SOURCE_DIR}/examples/step-*)

FOREACH(_step_full ${_steps})
  GET_FILENAME_COMPONENT(_step ${_step_full} NAME)

  FOREACH(_build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build_lowercase)

    SET(_test ${_category}/${_step}.${_build_lowercase})

    # Respect TEST_PICKUP_REGEX:
    IF( "${TEST_PICKUP_REGEX}" STREQUAL "" OR
        _test MATCHES "${TEST_PICKUP_REGEX}"  )

      SET(_step_dir ${CMAKE_CURRENT_BINARY_DIR}/${_step}.${_build_lowercase})

      FILE(GLOB _files ${_step_full}/*)
      SET(_command)
      FOREACH(_file ${_files})
        IF(NOT _file MATCHES "/(CMakeCache.txt|CMakeFiles|doc|.svn)")
          LIST(APPEND _command
            COMMAND ${CMAKE_COMMAND} -E copy ${_file} ${_step_dir}
            )
        ENDIF()
      ENDFOREACH()

      #
      # Add a "guard file" rule: The purpose of interrupt_guard is to force
      # a complete rerun of this test (CONFIGURE, BUILD and RUN stage) if
      # interrupt_guard is removed by run_test.cmake due to an
      # interruption.
      #
      ADD_CUSTOM_COMMAND(
        OUTPUT ${_step_dir}/interrupt_guard
        COMMAND ${CMAKE_COMMAND} -E make_directory ${_step_dir}
        COMMAND touch ${_step_dir}/interrupt_guard
        )

      # A rule how to copy the example step to the current directory:
      ADD_CUSTOM_COMMAND(OUTPUT ${_step_dir}/CMakeLists.txt
        ${_command}
        DEPENDS
          ${_files}
          ${DEAL_II_TARGET_DEBUG}
          ${DEAL_II_TARGET_RELEASE}
          ${_step_dir}/interrupt_guard
          )

      # And a rule on how to configure the example step:
      ADD_CUSTOM_COMMAND(OUTPUT ${_step_dir}/configure_output
        COMMAND rm -f ${_step_dir}/failing_configure_output
        COMMAND ${CMAKE_COMMAND}
          -DDEAL_II_DIR=${DEAL_II_BINARY_DIR} -DCMAKE_BUILD_TYPE=${_build} .
          > ${_step_dir}/configure_output 2>&1
          || (mv ${_step_dir}/configure_output
                 ${_step_dir}/failing_configure_output
              && echo "${_test}: CONFIGURE failed. Output:"
              && cat ${_step_dir}/failing_configure_output) # succeed anyway!
        WORKING_DIRECTORY ${_step_dir}
        DEPENDS ${_step_dir}/CMakeLists.txt
          ${DEAL_II_PATH}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PACKAGE_NAME}Config.cmake
        )

      # And a rule on how to build the example step:
      ADD_CUSTOM_COMMAND(OUTPUT ${_step_dir}/build_output
        COMMAND test ! -f ${_step_dir}/configure_output
          || (rm -f ${_step_dir}/failing_build_output
              &&${CMAKE_COMMAND} --build ${_step_dir} --target all
              > ${_step_dir}/build_output 2>&1)
          || (mv ${_step_dir}/build_output
                 ${_step_dir}/failing_build_output
              && echo "${_test}: CONFIGURE successful."
              && echo "${_test}: BUILD failed. Output:"
              && cat ${_step_dir}/failing_build_output
              && exit 1)
        COMMAND test -f ${_step_dir}/configure_output
          || (rm -f ${_step_dir}/build_output
              && rm -f ${_step_dir}/failing_build_output
              && echo "${_test}: BUILD stage not invoked due to failing CONFIGURE") # succeed anyway!
        WORKING_DIRECTORY ${_step_dir}
        DEPENDS ${_step_dir}/configure_output
        )

      # And a rule on how to run the example step:
      ADD_CUSTOM_COMMAND(OUTPUT ${_step_dir}/run_output
        COMMAND test ! -f ${_step_dir}/build_output
          || (rm -f ${_step_dir}/failing_run_output
             && ${CMAKE_COMMAND} --build ${_step_dir} --target run
             > ${_step_dir}/run_output 2>&1)
          || (mv ${_step_dir}/run_output
                 ${_step_dir}/failing_run_output
              && echo "${_test}: CONFIGURE successful."
              && echo "${_test}: BUILD successful."
              && echo "${_test}: RUN failed. Output:"
              && cat ${_step_dir}/failing_run_output
              && exit 1)
        COMMAND test -f ${_step_dir}/build_output
          || (   rm -f ${_step_dir}/run_output
              && rm -f ${_step_dir}/failing_run_output
              && echo "${_test}: RUN stage not invoked due to failing BUILD"
              && exit 1)
        WORKING_DIRECTORY ${_step_dir}
        DEPENDS ${_step_dir}/build_output
        )

      ITEM_MATCHES(_match ${_step} ${_${_build_lowercase}_steps})
      IF(_match)
        # Add a full test (CONFIGURE, BUILD, RUN):
        SET(_target ${_step}.${_build_lowercase}.run)
        ADD_CUSTOM_TARGET(${_target}
          COMMAND
               echo "${_test}: CONFIGURE successful."
            && echo "${_test}: BUILD successful."
            && echo "${_test}: RUN successful."
            && echo "${_test}: PASSED."
          DEPENDS ${_step_dir}/run_output
          )
      ELSE()
        # Add a minimal test (CONFIGURE, BUILD):
        SET(_target ${_step}.${_build_lowercase}.build)
        ADD_CUSTOM_TARGET(${_target}
          COMMAND
               echo "${_test}: CONFIGURE successful."
            && echo "${_test}: BUILD successful."
            && echo "${_test}: PASSED."
          DEPENDS ${_step_dir}/build_output
          )
      ENDIF()

      # And finally add the test:
      ADD_TEST(NAME ${_test}
        COMMAND ${CMAKE_COMMAND} -DTRGT=${_target} -DTEST=${_test}
          -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR}
          -DGUARD_FILE=${_step_dir}/interrupt_guard
          -P ${DEAL_II_SOURCE_DIR}/tests/run_test.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
      SET_TESTS_PROPERTIES(${_test} PROPERTIES
        LABEL "${_category}"
        TIMEOUT ${TEST_TIME_LIMIT}
        )

    ENDIF()
  ENDFOREACH()
ENDFOREACH()
