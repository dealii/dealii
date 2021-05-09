## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
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
# Tests whether it is possible to compile and link a dummy program with a
# given flag.
# If so, add it to variable.
#
# Usage:
#     ENABLE_IF_LINKS(variable flag)
#

MACRO(ENABLE_IF_LINKS _variable _flag)
  # keep on top to avoid cluttering the _flag and _flag_stripped variables
  ENABLE_IF_SUPPORTED(CMAKE_REQUIRED_FLAGS "-Werror")

  STRING(STRIP "${_flag}" _flag_stripped)
  IF(NOT "${_flag_stripped}" STREQUAL "")
    STRING(REGEX REPLACE "^-" "" _flag_name "${_flag_stripped}")
    STRING(REGEX REPLACE "\[-+,\]" "_" _flag_name "${_flag_name}")

    LIST(APPEND CMAKE_REQUIRED_LIBRARIES "${_flag_stripped}")
    CHECK_CXX_SOURCE_COMPILES("int main(){}" DEAL_II_HAVE_LINKER_FLAG_${_flag_name})
    RESET_CMAKE_REQUIRED()

    IF(DEAL_II_HAVE_LINKER_FLAG_${_flag_name})
      SET(${_variable} "${${_variable}} ${_flag_stripped}")
      STRING(STRIP "${${_variable}}" ${_variable})
    ENDIF()
  ENDIF()
ENDMACRO()
