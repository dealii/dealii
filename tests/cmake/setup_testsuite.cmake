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
# Setup testsuite:
#
# TODO: Describe and document the following:
#
# TEST_DIFF
# TEST_TIME_LIMIT
# NUMDIFF_DIR
#

#
# Wee need a diff tool, preferably numdiff:
#

FIND_PROGRAM(NUMDIFF_EXECUTABLE
  NAMES numdiff
  HINTS ${NUMDIFF_DIR}
  PATH_SUFFIXES bin
  )
FIND_PROGRAM(DIFF_EXECUTABLE
  NAMES diff
  )
IF( NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND"
    AND DIFF_EXECUTABLE MATCHES "-NOTFOUND" )
  MESSAGE(FATAL_ERROR
    "Could not find diff or numdiff. One of those are required for running the testsuite."
    )
ENDIF()

IF(NOT NUMDIFF_EXECUTABLE MATCHES "-NOTFOUND")
  SET_IF_EMPTY(TEST_DIFF "${NUMDIFF_EXECUTABLE} -a 1e-6 -q")
ELSE()
  SET_IF_EMPTY(TEST_DIFF "${DIFF_EXECUTABLE}")
ENDIF()
#
# Son, we have to talk about quotings:
#
SEPARATE_ARGUMENTS(TEST_DIFF UNIX_COMMAND ${TEST_DIFF})

MARK_AS_ADVANCED(DIFF_EXECUTABLE NUMDIFF_EXECUTABLE)

#
# Set a default time limit of 120 seconds:
#
SET_IF_EMPTY(TEST_TIME_LIMIT 120)

