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
# This is a small project that provides some custom targets that need a
# working, independent subproject because they are called back during
# configure phase.
#

IF(NOT DEFINED TEST_DIR)
  MESSAGE(FATAL_ERROR "\n
  The testsuite cannot be set up as an independent project.
  Please configure \"deal.II\" directly instead.\n\n"
    )
ENDIF()

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(testsuite NONE)

# Regenerate tests (run "make rebuild_cache" in subprojects):
ADD_CUSTOM_TARGET(regen_tests)

# Clean all tests
ADD_CUSTOM_TARGET(clean_tests)

FILE(GLOB _categories RELATIVE ${TEST_DIR} ${TEST_DIR}/*)
SET(_categories all-headers build_tests mesh_converter ${_categories})
LIST(REMOVE_DUPLICATES _categories)

FOREACH(_category ${_categories})
  IF( EXISTS ${CMAKE_SOURCE_DIR}/${_category}/CMakeLists.txt OR
      EXISTS ${TEST_DIR}/${_category}/CMakeLists.txt )

    ADD_CUSTOM_TARGET(regen_tests_${_category}
      COMMAND [ ! -d ${_category} ] || ${CMAKE_COMMAND}
        --build ${CMAKE_BINARY_DIR}/${_category} --target regenerate
      )
    ADD_DEPENDENCIES(regen_tests regen_tests_${_category})

    ADD_CUSTOM_TARGET(clean_tests_${_category}
      COMMAND [ ! -d ${_category} ] || ${CMAKE_COMMAND}
        --build ${CMAKE_BINARY_DIR}/${_category} --target clean
      )
    ADD_DEPENDENCIES(clean_tests clean_tests_${_category})

  ENDIF()
ENDFOREACH()
