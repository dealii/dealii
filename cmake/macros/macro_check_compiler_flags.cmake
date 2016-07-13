## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
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
# Usage:
#   CHECK_COMPILER_FLAGS(_compiler_flags_variable _linker_flags_variable _var)
#
# This macro tries to compile and link a simple "int main(){ return 0; }
# with the given set of compiler and flags provided in
# _compiler_flags_variable and _linker_flags_variable. If the test is
# succesful the variable ${_var} is set to true, otherwise to false.
#

MACRO(CHECK_COMPILER_FLAGS _compiler_flags_variable _linker_flags_variable _var)
  #
  # Rerun this test if flags have changed:
  #
  IF(NOT "${${_compiler_flags_variable}}" STREQUAL "${CACHED_${_var}_${_compiler_flags_variable}}"
     OR NOT "${${_linker_flags_variable}}" STREQUAL "${CACHED_${_var}_${_linker_flags_variable}}")
    UNSET(${_var} CACHE)
  ENDIF()

  SET(CACHED_${_var}_${_compiler_flags_variable} "${${_compiler_flags_variable}}"
    CACHE INTERNAL "" FORCE
    )
  SET(CACHED_${_var}_${_linker_flags_variable} "${${_linker_flags_variable}}"
    CACHE INTERNAL "" FORCE
    )

  SET(CMAKE_REQUIRED_FLAGS ${${_compiler_flags_variable}})
  SET(CMAKE_REQUIRED_LIBRARIES ${${_linker_flags_variable}})
  CHECK_CXX_SOURCE_COMPILES("int main(){ return 0; }" ${_var})
  RESET_CMAKE_REQUIRED()

  IF(${_var})
    SET(${_var} TRUE CACHE INTERNAL "")
  ENDIF()
ENDMACRO()
