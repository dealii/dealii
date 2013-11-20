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
# Setup the testsuite. The top level targets defined here merely act as a
# multiplexer for the ./tets/ project where all the actual work is
# done...
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
ADD_CUSTOM_TARGET(setup_tests
  COMMAND ${CMAKE_COMMAND}
    --build ${CMAKE_BINARY_DIR}/tests --target setup_tests
    -- ${MAKEOPTS}
  )

# Depend on a compiled library:
ADD_DEPENDENCIES(setup_tests build_library)

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

# Remove all tests:
ADD_CUSTOM_TARGET(prune_tests
  COMMAND ${CMAKE_COMMAND}
    --build ${CMAKE_BINARY_DIR}/tests --target prune_tests
    -- ${MAKEOPTS}
  )

#
# Setup the testsuite and pass all relevant "TEST_" and "_DIR" variables
# down to it:
#

MESSAGE(STATUS "Setup testsuite")

SET(_options)
LIST(APPEND _options -DDEAL_II_SOURCE_DIR=${CMAKE_SOURCE_DIR})
LIST(APPEND _options -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR})
FOREACH(_var
    DIFF_DIR NUMDIFF_DIR TEST_DIR TEST_DIFF TEST_OVERRIDE_LOCATION
    TEST_PICKUP_REGEX TEST_TIME_LIMIT
    )
  # always undefine:
  LIST(APPEND _options "-U${_var}")
  IF(DEFINED ${_var})
    LIST(APPEND _options "-D${_var}=${${_var}}")
  ENDIF()
ENDFOREACH()

EXECUTE_PROCESS(
  COMMAND ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${_options}
    ${CMAKE_SOURCE_DIR}/tests
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
  OUTPUT_QUIET
  )
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
