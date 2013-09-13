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
# Header tests are special:
#
# Construct a list of all header files and build up a test that just tries
# to compile a simple worker (test_header.cc) that only includes the given
# header file. We omit linking to safe some time.
#

SET(_category all-headers)

FILE(GLOB_RECURSE _header
  ${CMAKE_SOURCE_DIR}/include/*.h
  )

FOREACH(_full_file ${_header})

  GET_FILENAME_COMPONENT(_file ${_full_file} NAME)

  # TODO: A more sophisticated way to get the relative include path:
  GET_FILENAME_COMPONENT(_path ${_full_file} PATH)
  GET_FILENAME_COMPONENT(_path ${_path} NAME)
  IF("${_path}" STREQUAL "std_cxx1x")
    SET(_path "base/std_cxx1x")
  ENDIF()

  FOREACH(_build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build_lowercase)

    SET(_test ${_category}-${_path}-${_file}.${_build_lowercase})

    #
    # Respect TEST_PICKUP_REGEX:
    #
    IF( "${TEST_PICKUP_REGEX}" STREQUAL "" OR
        _test MATCHES "${TEST_PICKUP_REGEX}" )

      #
      # Add an object library for each header file and build configuration:
      #
      ADD_LIBRARY(${_test} OBJECT EXCLUDE_FROM_ALL test_header.cc)
      SET_TARGET_PROPERTIES(${_test} PROPERTIES
        COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
        COMPILE_FLAGS "${DEAL_II_CXX_FLAGS_${_build}}"
        )
      SET_PROPERTY(TARGET ${_test} APPEND PROPERTY
        INCLUDE_DIRECTORIES
          "${CMAKE_BINARY_DIR}/include"
          "${CMAKE_SOURCE_DIR}/include"
          "${CMAKE_SOURCE_DIR}/include/deal.II/"
        )
      SET_PROPERTY(TARGET ${_test} APPEND PROPERTY
        COMPILE_DEFINITIONS
          HEADER=<deal.II/${_path}/${_file}>
        )

      #
      # And finally add the test:
      #
      ADD_TEST(NAME ${_category}/${_path}/${_file}.${_build_lowercase}
        COMMAND ${CMAKE_COMMAND}
          -DTRGT=${_test}
          -DTEST=${_category}/${_path}/${_file}.${_build_lowercase}
          -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR}
          -P ${CMAKE_SOURCE_DIR}/cmake/scripts/run_test.cmake
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
      SET_TESTS_PROPERTIES(${_category}/${_path}/${_file}.${_build_lowercase} PROPERTIES
        LABEL "${_category}"
        TIMEOUT ${TEST_TIME_LIMIT}
        )
    ENDIF()
  ENDFOREACH()
ENDFOREACH()
