## ---------------------------------------------------------------------
##
## Copyright (C) 2016 - 2017 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Usage:
#   CHECK_COMPILER_SETUP("compiler flag string" "linker flag string" _var
#     [libraries]
#     )
#
# This macro tries to compile and link a simple "int main(){ return 0; }
# with the given set of provided compiler and linker flags  and an optional
# list of libraries to link against. If the test is successful the variable
# ${_var} is set to true, otherwise it is set to false.
#

MACRO(CHECK_COMPILER_SETUP _compiler_flags_unstr _linker_flags_unstr _var)
  #
  # Strip -Wl,--as-needed from the list of linker flags. This works around
  # a serious regression with ld.bfd in combination with -Wl,--as-needed
  # when compiling a simple
  #
  #    int main () { return 0; }
  #
  # and linking against a *huge* list of (entirely unused) libraries.
  #
  # Ideally, one should use ld.gold [1] instead of ld.bfd - but because
  # this is not always possible, simply disable -Wl,--as-needed.
  #
  # See https://github.com/dealii/dealii/issues/3686
  #
  # [1] https://lwn.net/Articles/274859/
  #
  STRING(REPLACE "-Wl,--as-needed" "" _linker_flags "${_linker_flags_unstr}")

  #
  # Strip leading and trailing whitespace to make CMake 2.8.8 happy
  #
  STRING(STRIP "${_compiler_flags_unstr}" _compiler_flags)
  STRING(STRIP "${_linker_flags}" _linker_flags)

  #
  # Rerun this test if flags have changed:
  #
  IF(NOT "${_compiler_flags}" STREQUAL "${CACHED_${_var}_compiler_flags}"
     OR NOT "${_linker_flags}" STREQUAL "${CACHED_${_var}_linker_flags}"
     OR NOT "${ARGN}" STREQUAL "${CACHED_${_var}_ARGN}")
    UNSET(${_var} CACHE)
  ENDIF()

  SET(CACHED_${_var}_compiler_flags "${_compiler_flags}"
    CACHE INTERNAL "" FORCE
    )
  SET(CACHED_${_var}_linker_flags "${_linker_flags}"
    CACHE INTERNAL "" FORCE
    )
  SET(CACHED_${_var}_ARGN "${ARGN}" CACHE INTERNAL "" FORCE)

  SET(CMAKE_REQUIRED_FLAGS ${_compiler_flags})
  SET(CMAKE_REQUIRED_LIBRARIES ${_linker_flags} ${ARGN})

  CHECK_CXX_SOURCE_COMPILES("int main(){ return 0; }" ${_var})
  RESET_CMAKE_REQUIRED()

  IF(${_var})
    SET(${_var} TRUE CACHE INTERNAL "")
  ENDIF()
ENDMACRO()
