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
# A Macro to pick up all tests in a test subdirectory
#
# If TEST_PICKUP_REGEX is set, only tests matching the regex will be
# processed.
#
# Usage:
#     DEAL_II_PICKUP_TESTS()
#

MACRO(DEAL_II_PICKUP_TESTS)
  GET_FILENAME_COMPONENT(_category ${CMAKE_CURRENT_SOURCE_DIR} NAME)

  FILE(GLOB _tests "*.output")

  FOREACH(_test ${_tests})

    #
    # Respect TEST_PICKUP_REGEX: Make sure we are allowed to pickup this
    # test:
    #
    IF( "${TEST_PICKUP_REGEX}" STREQUAL "" OR
        _test MATCHES "${TEST_PICKUP_REGEX}" )
      GET_FILENAME_COMPONENT(_test ${_test} NAME)
      SET(_define_test TRUE)
    ELSE()
      SET(_define_test FALSE)
    ENDIF()

    #
    # Query configuration and check whether we support it. Otherwise do
    # not define test:
    #
    STRING(REGEX MATCHALL "WITH_([0-9]|[A-Z]|_)*=(ON|OFF)" _matches ${_test})
    FOREACH(_match ${_matches})
      STRING(REGEX REPLACE "^(WITH_([0-9]|[A-Z]|_)*)=(ON|OFF)$" "\\1"
        _feature ${_match}
        )
      STRING(REGEX MATCH "(ON|OFF)$" _boolean ${_match})

      IF( (DEAL_II_${_feature} AND NOT ${_boolean}) OR
          (NOT DEAL_II_${_feature} AND ${_boolean}) )
        SET(_define_test FALSE)
      ENDIF()
    ENDFOREACH()

    IF(_define_test)
      #
      # TODO: mpirun support
      #

      IF(_test MATCHES debug)
        SET(_configuration DEBUG)
      ELSEIF(_test MATCHES release)
        SET(_configuration RELEASE)
      ELSE()
        SET(_configuration)
      ENDIF()

      SET(_comparison ${_test})
      STRING(REGEX REPLACE "\\..*" "" _test ${_test})
      DEAL_II_ADD_TEST(${_category} ${_test} ${_comparison} ${_configuration})
    ENDIF()

  ENDFOREACH()
ENDMACRO()

