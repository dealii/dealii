## ---------------------------------------------------------------------
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

MACRO(DEAL_II_PICKUP_TESTS)
  SET_IF_EMPTY(TEST_PICKUP_REGEX "$ENV{TEST_PICKUP_REGEX}")

  GET_FILENAME_COMPONENT(_category ${CMAKE_CURRENT_SOURCE_DIR} NAME)

  SET(DEAL_II_SOURCE_DIR) # avoid a bogus warning

  FILE(GLOB _tests "*.output")
  FOREACH(_test ${_tests})
    SET(_comparison ${_test})
    GET_FILENAME_COMPONENT(_test ${_test} NAME)

    #
    # Respect TEST_PICKUP_REGEX:
    #

    IF( "${TEST_PICKUP_REGEX}" STREQUAL "" OR
        "${_category}/${_test}" MATCHES "${TEST_PICKUP_REGEX}" )
      SET(_define_test TRUE)
    ELSE()
      SET(_define_test FALSE)
    ENDIF()

    #
    # Respect compiler constraint:
    #

    STRING(REGEX MATCHALL
      "compiler=[^=]*=(on|off|yes|no|true|false)" _matches ${_test}
      )
    FOREACH(_match ${_matches})
      STRING(REGEX REPLACE
        "^compiler=([^=]*)=(on|off|yes|no|true|false)$" "\\1"
        _compiler ${_match}
        )
      STRING(REGEX MATCH "(on|off|yes|no|true|false)$" _boolean ${_match})

      IF( ( "${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}"
              MATCHES "^${_compiler}"
            AND NOT ${_boolean} )
          OR ( NOT "${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}"
                   MATCHES "^${_compiler}"
               AND ${_boolean} ) )
        SET(_define_test FALSE)
      ENDIF()
    ENDFOREACH()

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

      # Make sure that _match is a valid feature constraint:
      IF(DEFINED DEAL_II_${_feature})
        STRING(REGEX MATCH "(on|off|yes|no|true|false)$" _boolean ${_match})
        IF( (DEAL_II_${_feature} AND NOT ${_boolean}) OR
            (NOT DEAL_II_${_feature} AND ${_boolean}) )
          SET(_define_test FALSE)
        ENDIF()
      ELSE()
        MESSAGE(FATAL_ERROR "
Invalid feature constraint \"${_match}\" in file
\"${_comparison}\":
The feature \"DEAL_II_${_feature}\" does not exist.\n"
          )
      ENDIF()
    ENDFOREACH()

    IF(_define_test)
      STRING(REGEX REPLACE "\\..*" "" _test ${_test})
      DEAL_II_ADD_TEST(${_category} ${_test} ${_comparison} ${_add_output})
    ENDIF()

  ENDFOREACH()
ENDMACRO()
