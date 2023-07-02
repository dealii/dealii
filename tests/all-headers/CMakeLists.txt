## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2023 by the deal.II authors
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

cmake_minimum_required(VERSION 3.13.4)
include(../scripts/setup_testsubproject.cmake)
project(testsuite CXX)

#
# Set up test environment (does not pick up anything):
#
deal_ii_pickup_tests()

#
# Header tests are special:
#
# Construct a list of all header files and build up a test that just tries
# to compile a simple worker (test_header.cc) that only includes the given
# header file. We omit linking to save some time.
#

set(_category all-headers)

#
# Glob together all header files and strip SOURCE_DIR/include/deal.II to
# get a correct relative path:
#
if(DEAL_II_BUILD_DIR)
  list(GET DEAL_II_INCLUDE_DIRS 1 _include_dir)
else()
  list(GET DEAL_II_INCLUDE_DIRS 0 _include_dir)
endif()
file(GLOB_RECURSE _headers RELATIVE ${_include_dir}/deal.II
  ${_include_dir}/deal.II/*.h
  )


# Do not test bundled headers to avoid issues when tests are run
# for an already installed library
string(REGEX REPLACE "bundled/[^;]+;?" "" _headers "${_headers}")

foreach(_header ${_headers})
  foreach(_build ${DEAL_II_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)

    set(_test ${_category}/${_header}.${_build_lowercase})
    string(REGEX REPLACE "\\/" "-" _target ${_header}.${_build_lowercase})

    # Respect TEST_PICKUP_REGEX:
    if( "${TEST_PICKUP_REGEX}" STREQUAL "" OR
        _test MATCHES "${TEST_PICKUP_REGEX}" )

      #
      # Add a "guard file" rule: The purpose of interrupt_guard.cc is to
      # force a complete rerun of this test (BUILD stage) if
      # interrupt_guard.cc is removed by run_test.cmake due to an
      # interruption.
      #
      add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_target}_interrupt_guard.cc
        COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/${_target}_interrupt_guard.cc
        )

      # Add an object library for each header file and build configuration:
      add_library(${_target} OBJECT EXCLUDE_FROM_ALL test_header.cc
        ${CMAKE_CURRENT_BINARY_DIR}/${_target}_interrupt_guard.cc
        )
      deal_ii_setup_target(${_target} ${_build})
      set_property(TARGET ${_target} APPEND PROPERTY
        COMPILE_DEFINITIONS HEADER=<deal.II/${_header}>
        )

      add_custom_target(${_target}.build
        COMMAND
             echo "${_test}: BUILD successful."
          && echo "${_test}: PASSED."
        )
      add_dependencies(${_target}.build ${_target})

      # And finally add the test:
      add_test(NAME ${_test}
        COMMAND ${CMAKE_COMMAND} -DTRGT=${_target}.build -DTEST=${_test}
          -DBINARY_DIR=${CMAKE_BINARY_DIR}
          -DGUARD_FILE=${CMAKE_CURRENT_BINARY_DIR}/${_target}_interrupt_guard.cc
          -P ${DEAL_II_PATH}/${DEAL_II_SHARE_RELDIR}/scripts/run_test.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
      set_tests_properties(${_test} PROPERTIES
        LABEL "${_category}"
        TIMEOUT ${TEST_TIME_LIMIT}
        )
    endif()
  endforeach()
endforeach()
