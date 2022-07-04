## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# This is a small worker to run a single test in the testsuite
#
# The following variables must be set:
#
#   TRGT - the name of the target that should be invoked
#   TEST - the test name (used for status messages)
#   BINARY_DIR - the build directory that contains the target
#
# Optional options:
#   EXPECT - the stage this test must reach to be considered successful
#            (return value 0)
#            Possible values are CONFIGURE, BUILD, RUN, DIFF, PASSED
#
#   GUARD_FILE - used to detect a forced interruption of this script: On
#                startup the backed up file ${GUARD_FILE}_bck is put back
#                in place as ${GUARD_FILE} and on exit ${GUARD_FILE} is
#                backed up as ${GUARD_FILE}_bck. If on startup a stale
#                ${GUARD_FILE} is found, it is deleted.
#

IF(NOT "${GUARD_FILE}" STREQUAL "" AND EXISTS ${GUARD_FILE})
  #
  # Guard file still exists, so this script must have been interrupted.
  # Remove guard file to force a complete rerun:
  #
  EXECUTE_PROCESS(COMMAND rm -f ${GUARD_FILE})
ELSEIF(NOT "${GUARD_FILE}" STREQUAL "" AND EXISTS ${GUARD_FILE}_bck)
  #
  # A backed up guard file exists. Put it back in place:
  #
  EXECUTE_PROCESS(COMMAND mv ${GUARD_FILE}_bck ${GUARD_FILE})
ENDIF()


IF("${EXPECT}" STREQUAL "")
  SET(EXPECT "PASSED")
ENDIF()

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
  --build . --target ${TRGT}
  WORKING_DIRECTORY ${BINARY_DIR}
  RESULT_VARIABLE _result_code # ignored ;-)
  OUTPUT_VARIABLE _output
  )

#
# Determine the stage a test reached: Possible values are
#   CONFIGURE  - the test started with a special configure stage and failed during configure
#   BUILD      - the test reached the build stage and a compilation error occured
#   RUN        - the test reached the run stage but the run terminated with an error
#   DIFF       - the test reached the diff stage but output differed
#   PASSED     - the test passed all stages
#

STRING(REGEX MATCH "${TEST}: CONFIGURE failed\\." _configure_regex "${_output}")
STRING(REGEX MATCH "${TEST}: BUILD failed\\." _build_regex "${_output}")
STRING(REGEX MATCH "${TEST}: RUN failed\\." _run_regex "${_output}")
STRING(REGEX MATCH "${TEST}: DIFF failed\\." _diff_regex "${_output}")
STRING(REGEX MATCH "${TEST}: PASSED\\." _passed_regex "${_output}")

IF(NOT "${_passed_regex}" STREQUAL "")
  SET(_stage PASSED)
ELSEIF(NOT "${_diff_regex}" STREQUAL "")
  SET(_stage DIFF)
ELSEIF(NOT "${_run_regex}" STREQUAL "")
  SET(_stage RUN)
ELSEIF(NOT "${_configure_regex}" STREQUAL "")
  SET(_stage CONFIGURE)
ELSE() # unconditionally, because "BUILD failed." doesn't have to be printed...
  SET(_stage BUILD)
ENDIF()

#
# Print out the test result:
#

MESSAGE("Test ${TEST}: ${_stage}")

MESSAGE("===============================   OUTPUT BEGIN  ===============================")

IF("${_stage}" STREQUAL "PASSED")
  STRING(REGEX REPLACE ".*\\/" "" _test ${TEST})
  #
  # MPI tests have a special runtime directory so rename:
  # test.mpirun=X.BUILD -> test.BUILD/mpirun=X
  #
  STRING(REGEX REPLACE "\\.(mpirun=[0-9]+)(\\..*)" "\\2/\\1" _test ${_test})
  #
  # Also output the diff file if we guessed the location correctly. This is
  # solely for cosmetic reasons: The diff file is either empty (if
  # comparison against the main comparison file was successful) or contains
  # a string explaining which comparison file variant succeeded.
  #
  SET(_diff "")
  IF(EXISTS ${BINARY_DIR}/${_test}/diff)
    FILE(READ ${BINARY_DIR}/${_test}/diff _diff)
  ENDIF()
  MESSAGE("${_diff}${TEST}: PASSED.")

ELSE()

  IF( "${_stage}" STREQUAL "BUILD" AND "${_build_regex}" STREQUAL "" )
    # Some special output in case the BUILD stage failed in a regression test:
    MESSAGE("${TEST}: BUILD failed. Output:")
  ENDIF()
  MESSAGE("${_output}")
  MESSAGE("")
  MESSAGE("${TEST}: ******    ${_stage} failed    *******")
  MESSAGE("")
ENDIF()

MESSAGE("===============================    OUTPUT END   ===============================")

#
# Back up guard file:
#

IF(NOT "${GUARD_FILE}" STREQUAL "" AND EXISTS ${GUARD_FILE})
  EXECUTE_PROCESS(COMMAND mv ${GUARD_FILE} ${GUARD_FILE}_bck)
ENDIF()

#
# Bail out:
#

IF(NOT "${_stage}" STREQUAL "${EXPECT}")
  MESSAGE("Expected stage ${EXPECT} - aborting")
  MESSAGE(FATAL_ERROR "*** abort")
ELSEIF(NOT "${_stage}" STREQUAL "PASSED")
  MESSAGE("Expected stage ${EXPECT} - test considered successful.")
ENDIF()
