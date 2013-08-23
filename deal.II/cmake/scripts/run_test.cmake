#
# TODO: The following variables have to be set:
#
# TEST
# DEAL_II_BINARY_DIR
#

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
  --build ${DEAL_II_BINARY_DIR} --target ${TEST}.diff
  RESULT_VARIABLE _result_code
  OUTPUT_VARIABLE _output
  )

IF("${_result_code}" STREQUAL "0")
  MESSAGE("***   ***")
  MESSAGE("${TEST}: Diff successful")
  MESSAGE(${_output})
  MESSAGE("***   ***")

ELSE()

  MESSAGE("***   ***")
  MESSAGE("${TEST}: Diff failed:")
  MESSAGE(${_output})
  MESSAGE("***   ***")
  MESSAGE(FATAL_ERROR "*** Test aborted.")
ENDIF()
