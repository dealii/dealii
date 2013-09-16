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

# Setup the testsuite:
INCLUDE(cmake/setup_testsuite.cmake)

#
# This is a bloody hack to avoid a severe performance regression with
# GNU Make that really does not like 12k top level targets at once...
#
# The only choice we have is to set up every test subdirectory as an
# independent project. Unfortunately this adds quite a significant amount
# of complexity :-(
#

MESSAGE(STATUS "Proceed to test definitions now.")

ADD_CUSTOM_TARGET(setup_test)
ADD_CUSTOM_TARGET(cleanup_test)
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/deal.IITestfile.cmake "")

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

FOREACH(_category ${_categories})
  # avoid stale test definitions:
  FILE(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/${_category}/CTestTestfile.cmake)

  ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/${_category}
    COMMAND ${CMAKE_COMMAND} -E make_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    )
  ADD_CUSTOM_TARGET(setup_test_${_category}
    COMMAND ${CMAKE_COMMAND} -E make_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    COMMAND ${CMAKE_COMMAND}
      -G${CMAKE_GENERATOR}
      -DDEAL_II_DIR=${CMAKE_BINARY_DIR}
      -DDEAL_II_SOURCE_DIR=${CMAKE_SOURCE_DIR}
      -DTEST_CMAKE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/cmake # TODO: Refactor
      -DTEST_DIFF="${TEST_DIFF}" # TODO: Refactor
      -DTEST_PICKUP_REGEX=${TEST_PICKUP_REGEX} # TODO: Refactor
      -DTEST_TIME_LIMIT=${TEST_TIME_LIMIT} # TODO: Refactor
      ${CMAKE_CURRENT_SOURCE_DIR}/${_category}
      > /dev/null # Shoo!
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_category}
    COMMENT "Processing ./tests/${_category}"
    )
  ADD_DEPENDENCIES(setup_test_${_category} library)
  ADD_DEPENDENCIES(setup_test setup_test_${_category})

  ADD_CUSTOM_TARGET(cleanup_test_${_category}
    COMMAND ${CMAKE_COMMAND} -E remove_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    )
  ADD_DEPENDENCIES(cleanup_test cleanup_test_${_category})

  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/deal.IITestfile.cmake
    "SUBDIRS(${_category})\n"
    )
ENDFOREACH()

SET_PROPERTY(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  APPEND PROPERTY TEST_INCLUDE_FILE deal.IITestfile.cmake
  )
