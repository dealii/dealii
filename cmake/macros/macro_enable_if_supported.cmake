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
# Tests whether the cxx compiler understands a flag.
# If so, add it to 'variable'.
#
# Usage:
#     enable_if_supported(variable flag)
#

macro(enable_if_supported _variable _flag)
  string(STRIP "${_flag}" _flag_stripped)

  #
  # Gcc does not emit a warning if testing -Wno-... flags which leads to
  # false positive detection. Unfortunately it later warns that an unknown
  # warning option is used if another warning is emitted in the same
  # compilation unit.
  # Therefore we invert the test for -Wno-... flags:
  #
  set(_flag_sanitized "${_flag_stripped}")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    string(REPLACE "-Wno-" "-W" _flag_sanitized "${_flag_stripped}")
  endif()

  if(NOT "${_flag_stripped}" STREQUAL "")
    string(REGEX REPLACE "^-" "" _flag_name "${_flag_stripped}")
    string(REGEX REPLACE "\[-+,\]" "_" _flag_name "${_flag_name}")

    CHECK_CXX_COMPILER_FLAG("${_flag_sanitized}" DEAL_II_HAVE_FLAG_${_flag_name})

    if(DEAL_II_HAVE_FLAG_${_flag_name})
      set(${_variable} "${${_variable}} ${_flag_stripped}")
      string(STRIP "${${_variable}}" ${_variable})
    endif()
  endif()
endmacro()

