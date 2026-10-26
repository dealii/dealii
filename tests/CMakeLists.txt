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

#
# Set up the testsuite.
#
# We define toplevel targets:
#    setup_tests    - set up testsuite subprojects
#    build_tests    - set up testsuite subprojects and build all test
#                     executables
#    prune_tests    - remove all testsuite subprojects
#    test           - run all quick_tests
#
# as well as setup_tests_<category>, build_tests_<category>, and
# prune_tests_<category> for every test category.
#

cmake_minimum_required(VERSION 3.13.4)

macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

set_if_empty(MAKEOPTS $ENV{MAKEOPTS})
set_if_empty(DEAL_II_DIR $ENV{DEAL_II_DIR})

#
# Setup for the case that tests are configured as a stand-alone project:
#

if(NOT DEFINED DEAL_II_HAVE_TESTS_DIRECTORY)
  message(STATUS "This is CMake ${CMAKE_VERSION}")
  message(STATUS "")

  if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "The testsuite cannot be configured in-source. "
      "Please create a separate build directory!"
      )
  endif()

  find_package(deal.II 9.8.0 REQUIRED HINTS ${DEAL_II_DIR})
  deal_ii_initialize_cached_variables()
  project(TESTSUITE CXX)
  set(_options "-DDEAL_II_DIR=${DEAL_II_PATH}")

  file(WRITE ${CMAKE_BINARY_DIR}/detailed.log
    "#        CMAKE_CXX_COMPILER:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}"
    )
endif()

#
# Setup for the case that tests are included via ADD_SUBDIRECTORY from the
# top-level CMake project:
#

if(DEFINED DEAL_II_HAVE_TESTS_DIRECTORY)
  message(STATUS "Setting up testsuite")

  #
  # Write minimalistic CTestTestfile.cmake files to CMAKE_BINARY_DIR:
  #
  file(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "subdirs(tests)\n")

  set(_options "-DDEAL_II_DIR=${CMAKE_BINARY_DIR}")
endif()

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake "")

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
# Configure the testsuite subprojects with the same generator as the main
# project:
#
list(APPEND _options "-G${CMAKE_GENERATOR}")
if(NOT "${CMAKE_GENERATOR_PLATFORM}" STREQUAL "")
  list(APPEND _options "-A${CMAKE_GENERATOR_PLATFORM}")
endif()
if(NOT "${CMAKE_GENERATOR_TOOLSET}" STREQUAL "")
  list(APPEND _options "-T${CMAKE_GENERATOR_TOOLSET}")
endif()
if(CMAKE_GENERATOR MATCHES "Makefiles|Ninja")
  list(APPEND _options "-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}")
endif()

#
# Find all testsuite subprojects, i.e., every directory that contains a
# CMakeLists.txt file.
#
set(_categories)
file(GLOB _dirs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*
  )
foreach(_dir ${_dirs})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_dir}/CMakeLists.txt)
    list(APPEND _categories ${_dir})
  endif()
endforeach()

#
# Custom targets for the testsuite:
#

#
# Clear the test summary files that records how many tests have been configured.
#
set(_summary_files_prefix "${CMAKE_CURRENT_BINARY_DIR}/test_summary")

# Set up all tests and print a summary:
add_custom_target(setup_tests
  COMMAND ${CMAKE_COMMAND}
    -D_summary_files_prefix="${_summary_files_prefix}"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/print_summary_files.cmake
  COMMENT "Setting up tests and printing test summary"
  )

# Set up all tests and build all test executables:
add_custom_target(build_tests)

# Remove all tests and also remove all test summary files:
add_custom_target(prune_tests
  COMMAND ${CMAKE_COMMAND}
    -D_summary_files_prefix="${_summary_files_prefix}"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/scripts/remove_summary_files.cmake
  COMMENT "Cleaning tests and removing test summary files"
  )

foreach(_category ${_categories})
  set(_category_dir ${CMAKE_CURRENT_SOURCE_DIR}/${_category})

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_category})

  if(DEAL_II_MSVC)
    set(_command ${CMAKE_COMMAND} ${_options} ${_category_dir})
  else()
    set(_command ${CMAKE_COMMAND} ${_options} ${_category_dir} > "${_summary_files_prefix}_${_category}")
  endif()

  add_custom_target(setup_tests_${_category}
    COMMAND ${_command}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    COMMENT "Processing tests/${_category}"
    )
  add_dependencies(setup_tests setup_tests_${_category})

  #
  # Build all test executables of the category. Test executables that fail
  # to build do not stop the build, their failure is reported by the
  # corresponding test. We run only one such build at a time (via the
  # console pool) because every build is already parallel.
  #
  add_custom_target(build_tests_${_category}
    COMMAND ${CMAKE_COMMAND}
      -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/${_category}
      -DGENERATOR=${CMAKE_GENERATOR}
      -DCONFIG=$<CONFIG>
      -P ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/scripts/build_tests.cmake
    COMMENT "Building tests/${_category}"
    USES_TERMINAL
    )
  add_dependencies(build_tests_${_category} setup_tests_${_category})
  add_dependencies(build_tests build_tests_${_category})

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

#
# Define a top-level "test" target that runs our quick tests wrapper.
#

add_custom_target(test
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/run_quick_tests.cmake
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Running quicktests..."
  USES_TERMINAL
  )

#
# The run_quick_tests.cmake script runs all tests in the a-framework,
# examples, and quick_tests categories. Make sure that they are set up and
# built:
#
foreach(_category a-framework examples quick_tests)
  if(TARGET build_tests_${_category})
    add_dependencies(test build_tests_${_category})
  endif()
endforeach()

if(DEFINED DEAL_II_HAVE_TESTS_DIRECTORY)
  #
  # Building the test executables requires the library targets and the
  # test driver:
  #
  foreach(_category ${_categories})
    foreach(_build ${DEAL_II_BUILD_TYPES})
      string(TOLOWER ${_build} _build_lowercase)
      add_dependencies(build_tests_${_category}
        ${DEAL_II_TARGET_NAME}_${_build_lowercase}
        )
    endforeach()
    add_dependencies(build_tests_${_category} deal_ii_test_driver)
  endforeach()

  #
  # Add a dummy target to make files known to IDEs like qtcreator
  #
  file(GLOB _misc
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/*/*.h
    )
  add_custom_target(testsuite SOURCES ${_misc})

  message(STATUS "Setting up testsuite - Done")
endif()
