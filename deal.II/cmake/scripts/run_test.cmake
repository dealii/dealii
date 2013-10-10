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
# This is a small worker to run a single test in the testsuite
#
# The following variables have to be set:
#
#   TRGT - the name of the target that should be invoked
#   TEST - the test name (used for status messages)
#   DEAL_II_BINARY_DIR - the build directory that contains the target
#

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
  --build ${DEAL_II_BINARY_DIR} --target ${TRGT}
  RESULT_VARIABLE _result_code
  OUTPUT_VARIABLE _output
  )

IF("${_result_code}" STREQUAL "0")

  MESSAGE("Test ${TEST}: PASSED")
  MESSAGE("===============================   OUTPUT BEGIN  ===============================")
  # Do not output everything, just that we are successful:
  MESSAGE("${TEST}: PASSED.")
  MESSAGE("===============================    OUTPUT END   ===============================")

ELSE()

  #
  # Determine whether the CONFIGURE, BUILD or RUN stages were run successfully:
  #

  # CONFIGURE is special because it only exists in build tests:
  STRING(REGEX MATCH "${TEST}: CONFIGURE successful\\." _configure_regex ${_output})
  STRING(REGEX MATCH "${TEST}: CONFIGURE failed\\." _configure_regex_fail ${_output})
  STRING(REGEX MATCH "${TEST}: BUILD successful\\." _build_regex ${_output})
  STRING(REGEX MATCH "${TEST}: RUN successful\\." _run_regex ${_output})
  IF(NOT "${_configure_regex_fail}" STREQUAL "")
    SET(_stage CONFIGURE)
  ELSEIF("${_build_regex}" STREQUAL "")
    SET(_stage BUILD)
  ELSEIF("${_run_regex}" STREQUAL "")
    SET(_stage RUN)
  ELSE()
    SET(_stage DIFF)
  ENDIF()

  MESSAGE("Test ${TEST}: ${_stage}")
  MESSAGE("===============================   OUTPUT BEGIN  ===============================")
  IF( "${_build_regex}" STREQUAL "" AND
      "${_configure_regex}" STREQUAL "" )
    # Some special output in case the BUILD stage failed in a regression test:
    MESSAGE("${TEST}: BUILD failed. Output:")
  ENDIF()
  MESSAGE(${_output})
  MESSAGE("")
  MESSAGE("${TEST}: ******    ${_stage} failed    *******")
  MESSAGE("")
  MESSAGE("===============================    OUTPUT END   ===============================")
  MESSAGE(FATAL_ERROR "*** abort")

ENDIF()
