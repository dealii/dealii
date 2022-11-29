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
# Set up the testsuite.
#
# We define toplevel targets:
#    setup_tests    - set up testsuite subprojects
#    prune_tests    - remove all testsuite subprojects
#

cmake_minimum_required(VERSION 3.3.0)

macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

set_if_empty(MAKEOPTS $ENV{MAKEOPTS})
set_if_empty(DEAL_II_DIR $ENV{DEAL_II_DIR})

#
# A bit of necessary setup if tests are configured as a stand-alone
# project:
#

if(NOT DEFINED DEAL_II_HAVE_TESTS_DIRECTORY)
  find_package(deal.II 9.3.0 REQUIRED HINTS ${DEAL_II_DIR})
  deal_ii_initialize_cached_variables()
  project(TESTSUITE CXX)
  file(WRITE ${CMAKE_BINARY_DIR}/detailed.log
    "#        CMAKE_CXX_COMPILER:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}"
    )
endif()

#
# A bit of necessary setup if tests are included via ADD_SUBDIRECTORY from
# the top-level CMake project:
#

if(DEFINED DEAL_II_HAVE_TESTS_DIRECTORY)

  #
  # We have to repeat the policy statement here because the new
  # CMAKE_MINIMUM_REQUIRED call resets our previous policy set in the main
  # CMakeLists.txt file.
  #
  if("${CMAKE_VERSION}" VERSION_LESS "3.11" AND POLICY CMP0037)
    # allow to override "test" target for quick tests
    cmake_policy(SET CMP0037 OLD)
  endif()

  #
  # If this CMakeLists.txt file is called from within the deal.II build
  # system, set up quick tests as well:
  #
  add_subdirectory(quick_tests)

  message(STATUS "Setting up testsuite")

  #
  # Write minimalistic CTestTestfile.cmake files to CMAKE_BINARY_DIR and
  # CMAKE_BINARY_DIR/tests:
  #
  file(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "subdirs(tests)\n")

  set(_options "-DDEAL_II_DIR=${CMAKE_BINARY_DIR}")

else()

  message(STATUS "This is CMake ${CMAKE_VERSION}")
  message(STATUS "")

  if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "The testsuite cannot be configured in-source. "
      "Please create a separate build directory!"
      )
  endif()

  find_package(deal.II 9.3.0 REQUIRED HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR})
  project(testsuite NONE)
  set(_options "-DDEAL_II_DIR=${DEAL_II_PATH}")
endif()

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
    "subdirs(quick_tests)\n"
    )

#
# Always undefine the following variables in the setup_tests target:
#
foreach(_var
    DIFF_DIR NUMDIFF_DIR TEST_PICKUP_REGEX TEST_TIME_LIMIT
    TEST_MPI_RANK_LIMIT TEST_THREAD_LIMIT ENABLE_PERFORMANCE_TESTS
    TESTING_ENVIRONMENT
    )
  list(APPEND _options "-U${_var}")
  if(NOT "${${_var}}" STREQUAL "")
    list(APPEND _options "-D${_var}=${${_var}}")
  endif()
endforeach()

#
# Find all testsuite subprojects, i.e., every directory that contains a
# CMakeLists.txt file (with the exception of "quick_tests").
#
set(_categories)
file(GLOB _dirs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*
  )
foreach(_dir ${_dirs})
  if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_dir}/CMakeLists.txt AND
      NOT ${_dir} MATCHES quick_tests)
    list(APPEND _categories ${_dir})
  endif()
endforeach()

#
# Custom targets for the testsuite:
#

# Setup tests:
add_custom_target(setup_tests)

# Remove all tests:
add_custom_target(prune_tests)

foreach(_category ${_categories})
  set(_category_dir ${CMAKE_CURRENT_SOURCE_DIR}/${_category})

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_category})

  if(DEAL_II_MSVC)
    set(_command ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${_options} ${_category_dir})
  else()
    # Do not pass the generator with -G so that we use make instead of ninja
    # for the test projects. This is because calling ninja several times in
    # parallel for the same project will break the configuration.
    set(_command ${CMAKE_COMMAND} ${_options} ${_category_dir} > /dev/null)
  endif()

  add_custom_target(setup_tests_${_category}
    COMMAND ${_command}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    COMMENT "Processing tests/${_category}"
    )
  add_dependencies(setup_tests setup_tests_${_category})

  add_custom_target(prune_tests_${_category}
    COMMAND ${CMAKE_COMMAND} -E remove_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    COMMAND ${CMAKE_COMMAND} -E make_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    COMMENT "Processing tests/${_category}"
    )
  add_dependencies(prune_tests prune_tests_${_category})

  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
    "subdirs(${_category})\n"
    )
endforeach()


if(DEFINED DEAL_II_HAVE_TESTS_DIRECTORY)
  #
  # Add a dummy target to make files known to IDEs like qtcreator
  #

  file(GLOB _misc
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.h
  )

  add_custom_target(dummy_tests_files
    SOURCES ${_misc}
  )

  message(STATUS "Setting up testsuite - Done")
endif()
