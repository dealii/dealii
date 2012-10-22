#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# TODO: Description...
#
MACRO(GET_CXX_SOURCE_RETURN_VALUE _source _var _exit_code)

  #
  # TODO: This file is still very basic :-]
  #

  IF(NOT DEFINED ${_var})
    FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx"
      "${_source}\n")

    MESSAGE(STATUS "Performing Test ${_var}")

    TRY_RUN(
      ${_var}_EXIT_CODE
      ${_var}_COMPILE_OK
      ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx)

    IF(${_var}_COMPILE_OK)
      SET(${_var} 1 CACHE INTERNAL "Test ${_var}")
      SET(${_exit_code} ${${_var}_EXIT_CODE} CACHE INTERNAL "Test ${_exit_code}")
      MESSAGE(STATUS "Performing Test ${_var} - Success")

      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
        "Performing C++ SOURCE FILE Test ${_var} succeded with the following output:\n"
        "${OUTPUT}\n"
        "Return value: ${${_var}}\n"
        "Source file was:\n${_source}\n")

    ELSE()

      SET(${_var} 0 CACHE INTERNAL "Test ${_var}")
      MESSAGE(STATUS "Performing Test ${_var} - Failed")

      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
        "Performing C++ SOURCE FILE Test ${_var} failed with the following output:\n"
        "${OUTPUT}\n"
        "Return value: ${${_var}_EXITCODE}\n"
        "Source file was:\n${_source}\n")
    ENDIF()
  ENDIF()
ENDMACRO()
