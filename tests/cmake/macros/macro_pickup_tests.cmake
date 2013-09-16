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
# A macro to pick up all tests in a test subdirectory
#
# If TEST_PICKUP_REGEX is set, only tests matching the regex will be
# processed.
#
# Usage:
#     DEAL_II_PICKUP_TESTS()
#

IF(NOT DEFINED SET_IF_EMPTY)
  MACRO(SET_IF_EMPTY _variable)
    IF("${${_variable}}" STREQUAL "")
      SET(${_variable} ${ARGN})
    ENDIF()
  ENDMACRO()
ENDIF()

MACRO(DEAL_II_PICKUP_TESTS)
  SET_IF_EMPTY(TEST_PICKUP_REGEX "$ENV{TEST_PICKUP_REGEX}")

  GET_FILENAME_COMPONENT(_category ${CMAKE_CURRENT_SOURCE_DIR} NAME)

  #
  # Strip the quoting from TEST_DIFF:
  #
  SET(_test_diff ${TEST_DIFF})
  SEPARATE_ARGUMENTS(TEST_DIFF UNIX_COMMAND ${TEST_DIFF})

  SET(DEAL_II_SOURCE_DIR) # avoid a bogus warning

  SET(_count "0")
  FILE(GLOB _tests "*.output")
  FOREACH(_test ${_tests})
    SET(_comparison ${_test})

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
    # Query configuration and check whether we support it. Otherwise
    # set _define_test to FALSE:
    #
    STRING(REGEX MATCHALL
      "with_([0-9]|[a-z]|_)*=(on|off|yes|no|true|false)" _matches ${_test}
      )
    FOREACH(_match ${_matches})
      STRING(REGEX REPLACE
        "^(with_([0-9]|[a-z]|_)*)=(on|off|yes|no|true|false)$" "\\1"
        _feature ${_match}
        )
      STRING(TOUPPER ${_feature} _feature)

      STRING(REGEX MATCH "(on|off|yes|no|true|false)$" _boolean ${_match})

      IF( (DEAL_II_${_feature} AND NOT ${_boolean}) OR
          (NOT DEAL_II_${_feature} AND ${_boolean}) )
        SET(_define_test FALSE)
      ENDIF()
    ENDFOREACH()

    IF(_define_test)
      #
      # Determine whether the test should be run with mpirun:
      #
      STRING(REGEX MATCH "mpirun=([0-9]*)" _n_cpu ${_test})
      IF("${_n_cpu}" STREQUAL "")
        # 0 indicates that no mpirun should be used
        SET(_n_cpu 0)
      ELSE()
        STRING(REGEX REPLACE "^mpirun=([0-9]*)$" "\\1" _n_cpu ${_n_cpu})
      ENDIF()

      IF(_test MATCHES debug)
        SET(_configuration DEBUG)
        MATH(EXPR _count "${_count}+1")
      ELSEIF(_test MATCHES release)
        SET(_configuration RELEASE)
        MATH(EXPR _count "${_count}+1")
      ELSE()
        SET(_configuration)
        MATH(EXPR _count "${_count}+2")
      ENDIF()

      STRING(REGEX REPLACE "\\..*" "" _test ${_test})
      DEAL_II_ADD_TEST(${_category} ${_test} ${_comparison} ${_n_cpu} ${_configuration})
    ENDIF()
  ENDFOREACH()

  MESSAGE(STATUS
    "Testsuite: Set up ${_count} regression tests in category ${_category}."
    )
  MESSAGE(STATUS
    "  (Timeout: ${TEST_TIME_LIMIT}, Diff: "${_test_diff}")"
    )
ENDMACRO()
