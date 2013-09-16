#
# TODO: The following variables have to be set:
#
# TARGET
# TEST
# DEAL_II_BINARY_DIR
#

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
  --build ${DEAL_II_BINARY_DIR} --target ${TRGT}
  RESULT_VARIABLE _result_code
  OUTPUT_VARIABLE _output
  )

IF("${_result_code}" STREQUAL "0")

  MESSAGE("Test ${TEST}: PASSED")
  MESSAGE("===============================   OUTPUT BEGIN  ===============================")
  #
  # Do not output everything, just that we are successful:
  #
  IF(TRGT MATCHES "\\.diff$") # ordinary test:
    MESSAGE("${TEST}: BUILD successful.")
    MESSAGE("${TEST}: RUN successful.")
    MESSAGE("${TEST}: DIFF successful.")
    MESSAGE("${TEST}: PASSED.")
  ELSEIF(TRGT MATCHES "\\.run$") # build_test for steps:
    MESSAGE("${TEST}: CONFIGURE successful.")
    MESSAGE("${TEST}: BUILD successful.")
    MESSAGE("${TEST}: RUN successful.")
    MESSAGE("${TEST}: PASSED.")
  ELSEIF(TRGT MATCHES "\\.build$") # build_test for steps:
    MESSAGE("${TEST}: CONFIGURE successful.")
    MESSAGE("${TEST}: BUILD successful.")
    MESSAGE("${TEST}: PASSED.")
  ELSE() # all-headers test:
    MESSAGE("${TEST}: BUILD successful.")
    MESSAGE("${TEST}: PASSED.")
  ENDIF()
  MESSAGE("===============================    OUTPUT END   ===============================")

ELSE()

  #
  # Determine whether the CONFIGURE, BUILD or RUN stages were run successfully:
  #

  # configure is special because it only exists in build tests:
  STRING(REGEX MATCH "${TEST}: CONFIGURE failed\\." _configure_regex ${_output})
  STRING(REGEX MATCH "${TEST}: BUILD successful\\." _build_regex ${_output})
  STRING(REGEX MATCH "${TEST}: RUN successful\\." _run_regex ${_output})
  IF(NOT "${_configure_regex}" STREQUAL "")
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
  IF("${_build_regex}" STREQUAL "")
    # Some special output in case the BUILD stage failed:
    MESSAGE("${TEST}: BUILD failed. Output:")
  ENDIF()
  MESSAGE(${_output})
  MESSAGE("")
  MESSAGE("${TEST}: ******    ${_stage} failed    *******")
  MESSAGE("")
  MESSAGE("===============================    OUTPUT END   ===============================")
  MESSAGE(FATAL_ERROR "*** abort")

ENDIF()
