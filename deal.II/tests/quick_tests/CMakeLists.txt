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
# A minimalistic set of tests:
#

ENABLE_TESTING()

FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER ${_build} _build_lowercase)

  #
  # Test whether thread affinity is well behaved:
  #
  SET(_target affinity.${_build_lowercase})

  ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL affinity.cc)
  DEAL_II_INSOURCE_SETUP_TARGET(${_target} ${_build})

  ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_target}-OK
    COMMAND ${_target} > ${CMAKE_CURRENT_BINARY_DIR}/${_target}-OK 2>&1 || (rm ${_target}-OK && exit 1)
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  ADD_CUSTOM_TARGET(${_target}.run DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_target}-OK)

  ADD_TEST(NAME ${_target}
    COMMAND ${CMAKE_COMMAND} -DTRGT=${_target}.run -DTEST=${_target}
      -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR}
      -P ${CMAKE_SOURCE_DIR}/cmake/scripts/run_test.cmake
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  SET_TESTS_PROPERTIES(${_target} PROPERTIES LABEL "sanity checks")

ENDFOREACH()

#
# A custom test target:
#
ADD_CUSTOM_TARGET(test
  COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
