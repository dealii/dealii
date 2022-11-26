## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
# Check for a compiler bug.
#
# Usage:
#     check_cxx_compiler_bug(source var),
#
# where source is a snipped of source code and var is a variable that will
# be set to true if the source could not be compiled and linked successfully.
# (This just inverts the logic of CHECK_CXX_SOURCE_COMPILES.)
#

macro(check_cxx_compiler_bug _source _var)
  if(NOT DEFINED ${_var}_OK)
    CHECK_CXX_SOURCE_COMPILES(
      "${_source}"
      ${_var}_OK
      )
    if(${_var}_OK)
      message(STATUS "Test successful, do not define ${_var}")
    else()
      message(STATUS "Test unsuccessful, define ${_var}")
    endif()
  endif()

  if(${_var}_OK)
    set(${_var})
  else()
    set(${_var} TRUE)
  endif()
endmacro()

