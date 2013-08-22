#
# TODO: The following variables have to be set:
#
# TEST
# TEST_SUFFIX
# DEAL_II_BINARY_DIR
#

MACRO(CALLBACK _target _target_suffix _msg_success _msg_error)
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
    --build ${DEAL_II_BINARY_DIR} --target ${_target}${_target_suffix}
    RESULT_VARIABLE _result_code
    OUTPUT_VARIABLE _output
    )

  IF("${_result_code}" STREQUAL "0")
    MESSAGE("${_target}: ${_msg_success}.")

  ELSE()

    MESSAGE("*** ${_msg_error}: ***")
    MESSAGE(${_output})
    MESSAGE(FATAL_ERROR "*** Test aborted.")
  ENDIF()
ENDMACRO()


CALLBACK(${TEST} "${TEST_SUFFIX}"
  "Build successful" "BUILD FAILED"
  )

CALLBACK(${TEST} ".run${TEST_SUFFIX}"
  "Run successful" "RUN FAILED"
  )

CALLBACK(${TEST} ".diff${TEST_SUFFIX}"
  "Diff successful" "DIFF FAILED"
  )

