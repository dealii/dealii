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
# Tests for the mesh_converter executable
#

SET(_category mesh_converter)

FILE(GLOB _meshes ${CMAKE_CURRENT_SOURCE_DIR}/meshes/?d/*.inp)

FOREACH(_full_file ${_meshes})
  GET_FILENAME_COMPONENT(_mesh ${_full_file} NAME_WE)

  SET(_test_full ${_category}/${_mesh})
  SET(_test_directory ${CMAKE_CURRENT_BINARY_DIR}/${_mesh}) # directory to run the test in
  STRING(REGEX REPLACE "\\.inp$" ".ucd" _comparison_file ${_full_file})

  # Is this a 3D or 2D mesh?
  SET(_dim 3)
  IF(_full_file MATCHES "meshes/2d")
    SET(_dim 2)
  ENDIF()

  # Respect TEST_PICKUP_REGEX:
  IF( "${TEST_PICKUP_REGEX}" STREQUAL "" OR
      _test_full MATCHES "${TEST_PICKUP_REGEX}" )

    FILE(MAKE_DIRECTORY ${_test_directory})

    ADD_CUSTOM_COMMAND(OUTPUT ${_test_directory}/output
      COMMAND rm -f ${_test_directory}/failing_output
      COMMAND touch ${_test_directory}/output
      COMMAND
        mesh_converter_exe ${_dim} ${_full_file} output
        || (mv ${_test_directory}/output
               ${_test_directory}/failing_output
            && echo "${_test_full}: RUN failed. ------ Result: ${_test_directory}/failing_output"
            && echo "${_test_full}: RUN failed. ------ Partial output:"
            && cat ${_test_directory}/failing_output
            && exit 1)
      COMMAND
        ${PERL_EXECUTABLE} -pi ${DEAL_II_SOURCE_DIR}/tests/normalize.pl
                               ${_test_directory}/output
      WORKING_DIRECTORY ${_test_directory}
      DEPENDS mesh_converter_exe ${DEAL_II_SOURCE_DIR}/tests/normalize.pl
      )
    ADD_CUSTOM_COMMAND(OUTPUT ${_test_directory}/diff
      COMMAND rm -f ${_test_directory}/failing_diff
      COMMAND touch ${_test_directory}/diff
      COMMAND
        ${TEST_DIFF} ${_comparison_file} ${_test_directory}/output > ${_test_directory}/diff
        || (mv ${_test_directory}/diff
               ${_test_directory}/failing_diff
            && echo "${_test_full}: RUN successful."
            && echo "${_test_full}: DIFF failed. ------ Source: ${_comparison_file}"
            && echo "${_test_full}: DIFF failed. ------ Result: ${_test_directory}/output"
            && echo "${_test_full}: DIFF failed. ------ Diff:   ${_test_directory}/failing_diff"
            && echo "${_test_full}: DIFF failed. ------ Diffs as follows:"
            && cat ${_test_directory}/failing_diff
            && exit 1)
      WORKING_DIRECTORY ${_test_directory}
      DEPENDS
        ${_test_directory}/output
        ${_comparison_file}
      )

    ADD_CUSTOM_TARGET(${_mesh}.diff DEPENDS ${_test_directory}/diff
      COMMAND
           echo "${_test_full}: RUN successful."
        && echo "${_test_full}: DIFF successful."
        && echo "${_test_full}: PASSED."
      )

    ADD_TEST(NAME ${_test_full}
      COMMAND ${CMAKE_COMMAND}
        -DTRGT=${_mesh}.diff
        -DTEST=${_test_full}
        -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR}
        -P ${DEAL_II_SOURCE_DIR}/tests/run_test.cmake
      WORKING_DIRECTORY ${_test_directory}
      )
    SET_TESTS_PROPERTIES(${_test_full} PROPERTIES
      LABEL "${_category}"
      TIMEOUT ${TEST_TIME_LIMIT}
      )

  ENDIF()

ENDFOREACH()
