## ---------------------------------------------------------------------
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
# Test whether a usable language compiler is available and if yes, call
# ENABLE_LANGUAGE(language)
#
# This works around a severe bug [1] in
#
#   ENABLE_LANGUAGE(Fortran OPTIONAL)
#
# [1] http://public.kitware.com/Bug/view.php?id=9220
#
# Usage:
#     ENABLE_LANGUAGE_FORTRAN_OPTIONAL(language)
#
# where language is either C or Fortran
#

MACRO(ENABLE_LANGUAGE_OPTIONAL _language)
  IF(NOT ${_language}_CHECKED)
    #
    # Run this check exactly once:
    #
    SET(${_language}_CHECKED TRUE CACHE INTERNAL "" FORCE)

    SET(_tmp ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/${_language}_test)
    file(REMOVE ${_tmp})

    IF(DEFINED CMAKE_${_language}_COMPILER)
      SET(_hint "-DCMAKE_${_language}_COMPILER=${CMAKE_${_language}_COMPILER}")
    ENDIF()

    FILE(WRITE ${_tmp}/CMakeLists.txt
      "PROJECT(foobar ${_language})"
      )

    IF(NOT "${CMAKE_TOOLCHAIN_FILE}" STREQUAL "")
      LIST(APPEND _hint "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
    ENDIF()

    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${_hint} .
      WORKING_DIRECTORY ${_tmp}
      RESULT_VARIABLE _result
      OUTPUT_QUIET
      ERROR_QUIET
      )

    IF("${_result}" STREQUAL "0")
      SET(DEAL_II_${_language}_COMPILER_WORKS TRUE CACHE INTERNAL "" FORCE)
      ENABLE_LANGUAGE(${_language})
    ELSE()
      MESSAGE(STATUS "No working ${_language} compiler found, disabling ${_language}")
    ENDIF()
  ELSE()
    #
    # Enable the language depending on the cached result from a former run:
    #
    IF(DEAL_II_${_language}_COMPILER_WORKS)
      ENABLE_LANGUAGE(${_language})
    ENDIF()
  ENDIF()
ENDMACRO()
