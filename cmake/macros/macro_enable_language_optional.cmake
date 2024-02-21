## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2022 by the deal.II authors
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
# Test whether a usable language compiler is available and if yes, call
# enable_language(language)
#
# This works around a severe bug [1] in
#
#   enable_language(Fortran OPTIONAL)
#
# [1] http://public.kitware.com/Bug/view.php?id=9220
#
# Usage:
#     enable_language_optional(language)
#
# where language is either C or Fortran
#

macro(enable_language_optional _language)
  if(NOT ${_language}_CHECKED)
    #
    # Run this check exactly once:
    #
    set(${_language}_CHECKED TRUE CACHE INTERNAL "" FORCE)

    set(_tmp ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/${_language}_test)
    file(REMOVE ${_tmp})

    if(DEFINED CMAKE_${_language}_COMPILER)
      set(_hint "-DCMAKE_${_language}_COMPILER=${CMAKE_${_language}_COMPILER}")
    endif()

    file(WRITE ${_tmp}/CMakeLists.txt
      "project(foobar ${_language})"
      )

    if(NOT "${CMAKE_TOOLCHAIN_FILE}" STREQUAL "")
      list(APPEND _hint "-DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
    endif()

    execute_process(
      COMMAND ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${_hint} .
      WORKING_DIRECTORY ${_tmp}
      RESULT_VARIABLE _result
      OUTPUT_QUIET
      ERROR_QUIET
      )

    if("${_result}" STREQUAL "0")
      set(DEAL_II_${_language}_COMPILER_WORKS TRUE CACHE INTERNAL "" FORCE)
      enable_language(${_language})
    else()
      message(STATUS "No working ${_language} compiler found, disabling ${_language}")
      set(DEAL_II_${_language}_COMPILER_WORKS FALSE CACHE INTERNAL "" FORCE)
    endif()
  else()
    #
    # Enable the language depending on the cached result from a former run:
    #
    if(${DEAL_II_${_language}_COMPILER_WORKS})
      enable_language(${_language})
    endif()
  endif()
endmacro()
