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

IF(NOT DEFINED DEAL_II_SOURCE_DIR OR NOT DEFINED DEAL_II_BINARY_DIR)
  MESSAGE(FATAL_ERROR "\n
  The testsuite cannot be set up as an independent project.
  Please configure \"deal.II\" directly instead.\n\n"
    )
ENDIF()

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(testsuite NONE)

#
# This is a bloody hack to avoid a severe performance penalty when using
# 12k top level targets with GNU Make that really does not like that...
#
# The only choice we have is to set up every test subdirectory as an
# independent project. Unfortunately this adds quite a significant amount
# of complexity :-(

#
# Custom targets to set and clean up the testsuite:
#

# Setup tests:
ADD_CUSTOM_TARGET(setup_tests)

# Regenerate tests (run "make rebuild_cache" in subprojects):
ADD_CUSTOM_TARGET(regen_tests)

# Clean all tests
ADD_CUSTOM_TARGET(clean_tests)

# Remove all tests:
ADD_CUSTOM_TARGET(prune_tests)

#
# Be as quiet as possible:
#

IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
  SET(_shoo)
ELSE()
  SET(_shoo >/dev/null)
ENDIF()

#
# Pass down all relevant configuration:
#

SET(_options)
FOREACH(_var
    DEAL_II_SOURCE_DIR DEAL_II_BINARY_DIR DIFF_DIR NUMDIFF_DIR
    TEST_DIFF TEST_OVERRIDE_LOCATION TEST_PICKUP_REGEX TEST_TIME_LIMIT
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
  IF(EXISTS ${CMAKE_SOURCE_DIR}/${_category}/CMakeLists.txt)
    SET(_category_dir ${CMAKE_SOURCE_DIR}/${_category})
  ELSEIF(EXISTS ${TEST_DIR}/${_category}/CMakeLists.txt)
    SET(_category_dir ${TEST_DIR}/${_category})
  ELSE()
    SET(_category_dir)
  ENDIF()

  IF(NOT "${_category_dir}" STREQUAL "")

    ADD_CUSTOM_TARGET(setup_tests_${_category}
      COMMAND ${CMAKE_COMMAND} -E make_directory
        ${CMAKE_BINARY_DIR}/${_category}
      COMMAND cd ${CMAKE_BINARY_DIR}/${_category} &&
        ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${_options} ${_category_dir}
        ${_shoo}
      DEPENDS ${_category_dir}
      COMMENT "Processing tests/${_category}"
      )
    ADD_DEPENDENCIES(setup_tests setup_tests_${_category})

    ADD_CUSTOM_TARGET(regen_tests_${_category}
      COMMAND [ ! -d ${_category} ] || ${CMAKE_COMMAND}
        --build ${CMAKE_BINARY_DIR}/${_category} --target regenerate
        ${_shoo}
      )
    ADD_DEPENDENCIES(regen_tests regen_tests_${_category})

    ADD_CUSTOM_TARGET(clean_tests_${_category}
      COMMAND [ ! -d ${_category} ] || ${CMAKE_COMMAND}
        --build ${CMAKE_BINARY_DIR}/${_category} --target clean
        ${_shoo}
      )
    ADD_DEPENDENCIES(clean_tests clean_tests_${_category})

    ADD_CUSTOM_TARGET(prune_tests_${_category}
      COMMAND ${CMAKE_COMMAND} -E remove_directory
        ${CMAKE_BINARY_DIR}/${_category}
      )
    ADD_DEPENDENCIES(prune_tests prune_tests_${_category})

    FILE(APPEND ${CMAKE_BINARY_DIR}/CTestTestfile.cmake
      "SUBDIRS(${_category})\n"
      )

  ENDIF()
ENDFOREACH()

