## ---------------------------------------------------------------------
## $Id$
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

#
# This is a bloody hack to avoid a severe performance penalty when using
# 12k top level targets with GNU Make that really does not like that...
#
# The only choice we have is to set up every test subdirectory as an
# independent project. Unfortunately this adds quite a significant amount
# of complexity :-(
#

#
# Setup the testsuite.
#

SET_IF_EMPTY(MAKEOPTS $ENV{MAKEOPTS})

MESSAGE(STATUS "")
MESSAGE(STATUS "Setup testsuite with TEST_DIR ${TEST_DIR}")

ADD_SUBDIRECTORY(
  ${CMAKE_SOURCE_DIR}/tests/quick_tests
  ${CMAKE_BINARY_DIR}/tests/quick_tests
  )

#
# Write a minimalistic CTestTestfile.cmake file to CMAKE_BINARY_DIR and
# CMAKE_BINARY_DIR/tests:
#

FILE(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake
  "SUBDIRS(tests)"
  )

FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
FILE(WRITE ${CMAKE_BINARY_DIR}/tests/CTestTestfile.cmake "")

#
# Custom targets to set and clean up the testsuite:
#

# Setup tests:
ADD_CUSTOM_TARGET(setup_tests)

# Remove all tests:
ADD_CUSTOM_TARGET(prune_tests)

# Regenerate tests (run "make rebuild_cache" in subprojects):
ADD_CUSTOM_TARGET(regen_tests
  COMMAND ${CMAKE_COMMAND}
    --build ${CMAKE_BINARY_DIR}/tests --target regen_tests
    -- ${MAKEOPTS}
  )

# Clean all tests
ADD_CUSTOM_TARGET(clean_tests
  COMMAND ${CMAKE_COMMAND}
    --build ${CMAKE_BINARY_DIR}/tests --target clean_tests
    -- ${MAKEOPTS}
  )

MESSAGE(STATUS "Setup testsuite")

#
# Provide custom targets to setup and prune the testsuite subproject:
#

EXECUTE_PROCESS(
  COMMAND ${CMAKE_COMMAND} -G${CMAKE_GENERATOR}
    -DTEST_DIR=${TEST_DIR}
    ${CMAKE_SOURCE_DIR}/tests
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
  OUTPUT_QUIET
  )

SET(_options)
LIST(APPEND _options -DDEAL_II_SOURCE_DIR=${CMAKE_SOURCE_DIR})
LIST(APPEND _options -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR})
FOREACH(_var
    DIFF_DIR
    NUMDIFF_DIR
    TEST_DIFF
    TEST_PICKUP_REGEX
    TEST_TIME_LIMIT
    MPIEXEC
    MPIEXEC_NUMPROC_FLAG
    MPIEXEC_PREFLAGS
    MPIEXEC_POSTFLAGS
    )
  # always undefine:
  LIST(APPEND _options "-U${_var}")
  IF(DEFINED ${_var})
    LIST(APPEND _options "-D${_var}=${${_var}}")
  ENDIF()
ENDFOREACH()

#
# Glob together a list of all subfolders to set up:
#

FILE(GLOB _categories RELATIVE ${TEST_DIR} ${TEST_DIR}/*)
SET(_categories all-headers build_tests mesh_converter ${_categories})
LIST(REMOVE_DUPLICATES _categories)

#
# Define a subproject for every enabled category:
#

FOREACH(_category ${_categories})
  IF(EXISTS ${CMAKE_SOURCE_DIR}/tests/${_category}/CMakeLists.txt)
    SET(_category_dir ${CMAKE_SOURCE_DIR}/tests/${_category})
  ELSEIF(EXISTS ${TEST_DIR}/${_category}/CMakeLists.txt)
    SET(_category_dir ${TEST_DIR}/${_category})
  ELSE()
    SET(_category_dir)
  ENDIF()

  IF(NOT "${_category_dir}" STREQUAL "")
    ADD_CUSTOM_TARGET(setup_tests_${_category}
      COMMAND ${CMAKE_COMMAND} -E make_directory
        ${CMAKE_BINARY_DIR}/tests/${_category}
      COMMAND cd ${CMAKE_BINARY_DIR}/tests/${_category} &&
        ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${_options} ${_category_dir}
        > /dev/null
      DEPENDS ${_category_dir}
      COMMENT "Processing tests/${_category}"
      )
    ADD_DEPENDENCIES(setup_tests setup_tests_${_category})

    ADD_CUSTOM_TARGET(prune_tests_${_category}
      COMMAND ${CMAKE_COMMAND} -E remove_directory
        ${CMAKE_BINARY_DIR}/tests/${_category}
      )
    ADD_DEPENDENCIES(prune_tests prune_tests_${_category})

    FILE(APPEND ${CMAKE_BINARY_DIR}/tests/CTestTestfile.cmake
      "SUBDIRS(${_category})\n"
      )
  ENDIF()
ENDFOREACH()
MESSAGE(STATUS "Setup testsuite - Done")

MESSAGE(STATUS "Regenerating testsuite subprojects")
EXECUTE_PROCESS(
  COMMAND ${CMAKE_COMMAND}
    --build ${CMAKE_BINARY_DIR}/tests --target regen_tests
    -- ${MAKEOPTS}
  OUTPUT_QUIET
  )
MESSAGE(STATUS "Regenerating testsuite subprojects - Done")

MESSAGE(STATUS "")
