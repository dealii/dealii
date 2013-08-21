#
# The following variables have to be set:
#
# TEST
# DEAL_II_BINARY_DIR
#

MACRO(CALLBACK _target _msg_success _msg_error)
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
    --build ${DEAL_II_BINARY_DIR} --target ${_target}/fast
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


CALLBACK(${TEST}
  "Build successful" "BUILD FAILED"
  )

CALLBACK(${TEST}.run
  "Run successful" "RUN FAILED"
  )

CALLBACK(${TEST}.diff
  "Diff successful" "DIFF FAILED"
  )

