## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2013 - 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.13.4)
include(../scripts/setup_testsubproject.cmake)
project(testsuite CXX)

#
# Set up test environment (does not pick up anything):
#

deal_ii_pickup_tests()

#
# Include deal.II targets:
#

if(NOT DEAL_II_TARGET_CONFIG_INCLUDED)
  include(${DEAL_II_TARGET_CONFIG})
  set(DEAL_II_TARGET_CONFIG_INCLUDED TRUE)
endif()

#
# Header tests are special:
#
# Construct a list of all header files and build up a test that just tries
# to compile a simple worker (test_header.cc) that only includes the given
# header file. We omit linking to save some time.
#

if(NOT TARGET compile_test_executables)
  add_custom_target(compile_test_executables)
endif()

set(_category all-headers)

# Do not test bundled headers to avoid issues when tests are run
# for an already installed library
string(REGEX REPLACE "bundled/[^;]+;?" "" _headers "${_headers}")

foreach(_build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER ${_build} _build_lowercase)

  #
  # Extract all include directories from the imported deal.II target.
  #

  get_target_property(_include_dirs ${DEAL_II_TARGET_${_build}}
    INTERFACE_INCLUDE_DIRECTORIES
    )

  #
  # Glob together all header files and strip ".../include/deal.II" to get a
  # correct relative path. By convention the first (and in case of a binary
  # directory also second) directory is ours:
  #

  list(GET _include_dirs 0 _include_dir)
  file(GLOB_RECURSE _headers RELATIVE ${_include_dir}/deal.II
    ${_include_dir}/deal.II/*.h
    )

  if(DEAL_II_BUILD_DIR)
    list(GET _include_dirs 1 _include_dir)
    file(GLOB_RECURSE _source_headers RELATIVE ${_include_dir}/deal.II
      ${_include_dir}/deal.II/*.h
      )
    list(APPEND _headers ${_source_headers})
  endif()

  #
  # Set up a test for each header file:
  #

  foreach(_header ${_headers})

    set(_test ${_category}/${_header}.${_build_lowercase})
    string(REGEX REPLACE "\\/" "-" _target ${_header}.${_build_lowercase})

    # Respect TEST_PICKUP_REGEX:
    if( "${TEST_PICKUP_REGEX}" STREQUAL "" OR
        _test MATCHES "${TEST_PICKUP_REGEX}" )

      #
      # Add an object library for each header file and build
      # configuration. Compile commands are run through the test driver
      # that records the output of failing commands in a build log
      # directory (see deal_ii_add_test()).
      #
      set(_build_log_dir "${CMAKE_CURRENT_BINARY_DIR}/build-logs/${_target}")

      add_library(${_target} OBJECT EXCLUDE_FROM_ALL test_header.cc)
      deal_ii_setup_target(${_target} ${_build})
      set_property(TARGET ${_target} APPEND PROPERTY
        COMPILE_DEFINITIONS HEADER=<deal.II/${_header}>
        )
      set_target_properties(${_target} PROPERTIES
        RULE_LAUNCH_COMPILE
          "\"${DEAL_II_TEST_DRIVER}\" launch --log-dir \"${_build_log_dir}\" --name <OBJECT> --"
        )
      add_dependencies(compile_test_executables ${_target})

      # And finally add the test:
      add_test(NAME ${_test}
        COMMAND ${DEAL_II_TEST_DRIVER} run
          --test ${_test}
          --log-dir ${_build_log_dir}
          --require "$<TARGET_OBJECTS:${_target}>"
          --build-only
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
      file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/test_targets.txt
        "${_test} ${_target} .\n"
        )
      set_tests_properties(${_test} PROPERTIES
        LABEL "${_category}"
        TIMEOUT ${TEST_TIME_LIMIT}
        )
    endif()
  endforeach()
endforeach()
