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

  MESSAGE("${TEST}: Test successful:")
  MESSAGE("===============================   OUTPUT BEGIN  ===============================")
  # Do not output everything, just that we are successful:
  MESSAGE("${TEST}: BUILD successful")
  MESSAGE("${TEST}: RUN successful")
  MESSAGE("${TEST}: DIFF successful")
  MESSAGE("===============================    OUTPUT END   ===============================")

ELSE()

  MESSAGE("${TEST}: Test failed:")
  MESSAGE("===============================   OUTPUT BEGIN  ===============================")
  MESSAGE(${_output})
  MESSAGE("===============================    OUTPUT END   ===============================")
  MESSAGE(FATAL_ERROR "*** abort")

ENDIF()
