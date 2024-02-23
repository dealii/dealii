## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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

if("${EXPECT}" STREQUAL "")
  set(EXPECT "PASSED")
endif()

execute_process(COMMAND ${CMAKE_COMMAND}
  --build . --target ${TRGT}
  WORKING_DIRECTORY ${BINARY_DIR}
  RESULT_VARIABLE _result_code # ignored ;-)
  OUTPUT_VARIABLE _output
  )

#
# Determine the stage a test reached: Possible values are
#   CONFIGURE  - the test started with a special configure stage and failed during configure
#   BUILD      - the test reached the build stage and a compilation error occurred
#   RUN        - the test reached the run stage but the run terminated with an error
#   DIFF       - the test reached the diff stage but output differed
#   PASSED     - the test passed all stages
#

string(REGEX MATCH "${TEST}: CONFIGURE failed\\." _configure_regex "${_output}")
string(REGEX MATCH "${TEST}: BUILD failed\\." _build_regex "${_output}")
string(REGEX MATCH "${TEST}: RUN failed\\." _run_regex "${_output}")
string(REGEX MATCH "${TEST}: DIFF failed\\." _diff_regex "${_output}")
string(REGEX MATCH "${TEST}: PASSED\\." _passed_regex "${_output}")

if(NOT "${_passed_regex}" STREQUAL "")
  set(_stage PASSED)
elseif(NOT "${_diff_regex}" STREQUAL "")
  set(_stage DIFF)
elseif(NOT "${_run_regex}" STREQUAL "")
  set(_stage RUN)
elseif(NOT "${_configure_regex}" STREQUAL "")
  set(_stage CONFIGURE)
else() # unconditionally, because "BUILD failed." doesn't have to be printed...
  set(_stage BUILD)
endif()

#
# Print out the test result:
#

message("Test ${TEST}: ${_stage}")

message("===============================   OUTPUT BEGIN  ===============================")

if("${_stage}" STREQUAL "PASSED")
  string(REGEX REPLACE ".*\\/" "" _test ${TEST})
  #
  # MPI tests have a special runtime directory so rename:
  # test.mpirun=X.BUILD -> test.BUILD/mpirun=X
  #
  string(REGEX REPLACE "\\.(mpirun=[0-9]+)(\\..*)" "\\2/\\1" _test ${_test})
  #
  # Also output the diff file if we guessed the location correctly. This is
  # solely for cosmetic reasons: The diff file is either empty (if
  # comparison against the main comparison file was successful) or contains
  # a string explaining which comparison file variant succeeded.
  #
  set(_diff "")
  if(EXISTS ${BINARY_DIR}/${_test}/diff)
    file(READ ${BINARY_DIR}/${_test}/diff _diff)
  endif()
  message("${_diff}${TEST}: PASSED.")

else()

  if( "${_stage}" STREQUAL "BUILD" AND "${_build_regex}" STREQUAL "" )
    # Some special output in case the BUILD stage failed in a regression test:
    message("${TEST}: BUILD failed. Output:")
  endif()
  message("${_output}")
  message("")
  message("${TEST}: ******    ${_stage} failed    *******")
  message("")
endif()

message("===============================    OUTPUT END   ===============================")

#
# Bail out:
#

if(NOT "${_stage}" STREQUAL "${EXPECT}")
  message("Expected stage ${EXPECT} - aborting")
  message(FATAL_ERROR "*** abort")
elseif(NOT "${_stage}" STREQUAL "PASSED")
  message("Expected stage ${EXPECT} - test considered successful.")
endif()
