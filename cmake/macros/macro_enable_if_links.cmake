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
#     enable_if_links(variable flag)
#

macro(enable_if_links _variable _flag)
  # keep on top to avoid cluttering the _flag and _flag_stripped variables
  enable_if_supported(CMAKE_REQUIRED_FLAGS "-Werror")

  string(STRIP "${_flag}" _flag_stripped)
  if(NOT "${_flag_stripped}" STREQUAL "")
    string(REGEX REPLACE "^-" "" _flag_name "${_flag_stripped}")
    string(REGEX REPLACE "\[-+,\]" "_" _flag_name "${_flag_name}")

    list(APPEND CMAKE_REQUIRED_LIBRARIES "${_flag_stripped}")
    CHECK_CXX_SOURCE_COMPILES("int main(){}" DEAL_II_HAVE_LINKER_FLAG_${_flag_name})
    reset_cmake_required()

    if(DEAL_II_HAVE_LINKER_FLAG_${_flag_name})
      set(${_variable} "${${_variable}} ${_flag_stripped}")
      string(STRIP "${${_variable}}" ${_variable})
    endif()
  endif()
endmacro()
