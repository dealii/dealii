#
# TODO: Description...
#
MACRO(GET_CXX_SOURCE_RETURN_VALUE SOURCE VAR EXIT_CODE)

  #
  # TODO: This file is still very basic :-]
  #

  FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx"
    "${SOURCE}\n")

  MESSAGE(STATUS "Performing Test ${VAR}")

  TRY_RUN(
    ${VAR}_EXIT_CODE
    ${VAR}_COMPILE_OK
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx)

  IF(${VAR}_COMPILE_OK)
    SET(${VAR} 1 CACHE INTERNAL "Test ${VAR}")
    SET(${EXIT_CODE} ${VAR_EXIT_CODE} CACHE INTERNAL "Test ${EXIT_CODE}")
    MESSAGE(STATUS "Performing Test ${VAR} - Success")

    FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Performing C++ SOURCE FILE Test ${VAR} succeded with the following output:\n"
      "${OUTPUT}\n"
      "Return value: ${${VAR}}\n"
      "Source file was:\n${SOURCE}\n")

  ELSE()

    SET(${VAR} 0 CACHE INTERNAL "Test ${VAR}")
    MESSAGE(STATUS "Performing Test ${VAR} - Failed")

    FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Performing C++ SOURCE FILE Test ${VAR} failed with the following output:\n"
      "${OUTPUT}\n"
      "Return value: ${${VAR}_EXITCODE}\n"
      "Source file was:\n${SOURCE}\n")
  ENDIF()

ENDMACRO()
