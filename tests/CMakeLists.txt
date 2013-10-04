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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
MESSAGE(STATUS "Proceed to test definitions now.")
ENABLE_TESTING()

IF(NOT "${CMAKE_PROJECT_NAME}" STREQUAL "deal.II")
  MESSAGE(FATAL_ERROR
    "The testsuite cannot be set up as an independent project. Please configure \"deal.II\" instead."
    )
ENDIF()

#
# Write a minimalistic CTestTestfile.cmake file to CMAKE_BINARY_DIR:
#
FILE(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake
  "SUBDIRS(tests)"
  )

#
# The hook into CTest so that the tests from the suprojects get picked up
# by the main project:
#

FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/deal.IITestfile.cmake "")
SET_PROPERTY(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  APPEND PROPERTY TEST_INCLUDE_FILE deal.IITestfile.cmake
  )

#
# Pass all relevant "TEST_" and "_DIR" variables down to the subprojects:
#

SET(_options)

LIST(APPEND _options -DDEAL_II_SOURCE_DIR=${CMAKE_SOURCE_DIR})
LIST(APPEND _options -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR})

FOREACH(_var TEST_DIR TEST_DIFF TEST_TIME_LIMIT TEST_PICKUP_REGEX NUMDIFF_DIR)
  LIST(APPEND _options "-U${_var}")
  IF(DEFINED ${_var})
    LIST(APPEND _options "-D${_var}=${${_var}}")
  ENDIF()
ENDFOREACH()

#
# Custom targets to set and clean up the testsuite:
#

# Setup tests:
ADD_CUSTOM_TARGET(setup_test)

# Clean all tests
ADD_CUSTOM_TARGET(clean_test)

# Remove all tests:
ADD_CUSTOM_TARGET(prune_test)

#
# The list of subfolders to set up:
#

SET(_categories
  all-headers aniso base bits build_tests codim_one deal.II fe grid hp integrators lac matrix_free multigrid serialization
  )
IF(DEAL_II_WITH_LAPACK)
  LIST(APPEND _categories lapack)
ENDIF()
IF(DEAL_II_WITH_METIS)
  LIST(APPEND _categories metis)
ENDIF()
IF(DEAL_II_WITH_MPI AND DEAL_II_WITH_P4EST)
  LIST(APPEND _categories gla mpi)
ENDIF()
IF(DEAL_II_WITH_P4EST)
  LIST(APPEND _categories distributed_grids)
ENDIF()
IF(DEAL_II_WITH_PETSC)
  LIST(APPEND _categories petsc)
ENDIF()
IF(DEAL_II_WITH_SLEPC)
  LIST(APPEND _categories slepc)
ENDIF()
IF(DEAL_II_WITH_TRILINOS)
  LIST(APPEND _categories trilinos)
ENDIF()
IF(DEAL_II_WITH_UMFPACK)
  LIST(APPEND _categories umfpack)
ENDIF()

#
# Define a subproject for every enabled category:
#

FOREACH(_category ${_categories})
  ADD_CUSTOM_TARGET(setup_test_${_category}
    COMMAND ${CMAKE_COMMAND} -E make_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    COMMAND cd ${CMAKE_CURRENT_BINARY_DIR}/${_category} &&
      ${CMAKE_COMMAND} -G${CMAKE_GENERATOR}
      ${_options}
      ${CMAKE_CURRENT_SOURCE_DIR}/${_category}
      > /dev/null # Shoo!
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_category}
    COMMENT "Processing ./tests/${_category}"
    )
  ADD_DEPENDENCIES(setup_test setup_test_${_category})

  # depend on a valid build directory (libraries built, config in place):
  ADD_DEPENDENCIES(setup_test_${_category} setup_build_dir)

  ADD_CUSTOM_TARGET(clean_test_${_category}
    COMMAND [ ! -d ${_category} ] || ${CMAKE_COMMAND}
      --build ${CMAKE_CURRENT_BINARY_DIR}/${_category} --target clean
    )
  ADD_DEPENDENCIES(clean_test clean_test_${_category})

  ADD_CUSTOM_TARGET(prune_test_${_category}
    COMMAND ${CMAKE_COMMAND} -E remove_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    )
  ADD_DEPENDENCIES(prune_test prune_test_${_category})

  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/deal.IITestfile.cmake
    "SUBDIRS(${_category})\n"
    )
ENDFOREACH()
