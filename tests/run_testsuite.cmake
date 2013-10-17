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
# Sanity check: Do not run the testsuite within a test directory
#

IF( "${CTEST_BINARY_DIRECTORY}" STREQUAL ""
    AND "${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "${CMAKE_CURRENT_LIST_DIR}" )
  MESSAGE(FATAL_ERROR "
ctest was invoked in the source directory (or test source directory) and
CTEST_BINARY_DIRECTORY is not set. Please either call ctest from within a
designated build directory, or set CTEST_BINARY_DIRECTORY accordingly.
"
    )
ENDIF()

#
# Try to find the source directory and invoke
# ./cmake/scripts/run_testsuite.cmake from this location:
#

IF("${CTEST_SOURCE_DIRECTORY}" STREQUAL "")
  # If CTEST_SOURCE_DIRECTORY is not set we just assume that this script
  # was called residing under ../tests relative to the source directory.
  GET_FILENAME_COMPONENT(_path "${CMAKE_CURRENT_LIST_DIR}" PATH)
  SET(CTEST_SOURCE_DIRECTORY ${_path}/deal.II)
ENDIF()

IF(NOT EXISTS ${CTEST_SOURCE_DIRECTORY}/cmake/scripts/run_testsuite.cmake)
  MESSAGE(FATAL_ERROR "
Could not find a suitable source directory.
There is no source directory \"../deal.II\" relative to the location of
this script. Please, set CTEST_SOURCE_DIRECTORY manually to the appropriate
source directory.
"
    )
ENDIF()

MESSAGE("-- Redirect to: ${CTEST_SOURCE_DIRECTORY}/cmake/scripts/run_testsuite.cmake")
INCLUDE(${CTEST_SOURCE_DIRECTORY}/cmake/scripts/run_testsuite.cmake)
