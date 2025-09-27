## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
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
# Tests whether the cxx compiler understands a flag. If so, add it to
# 'variable'.
#
# Note: This macro will reset the CMAKE_REQUIRED_* variables.
#
# Usage:
#     enable_if_supported(variable flag)
#

macro(enable_if_supported _variable _flag)
  reset_cmake_required()

  #
  # Add -Werror if available:
  #
  CHECK_CXX_COMPILER_FLAG("-Werror" DEAL_II_HAVE_FLAG_Werror)
  if(DEAL_II_HAVE_FLAG_Werror)
    string(STRIP "${CMAKE_REQUIRED_FLAGS} -Werror" CMAKE_REQUIRED_FLAGS)
  endif()

  #
  # Gcc does not emit a warning if testing -Wno-... flags which leads to
  # false positive detection. Unfortunately it later warns that an unknown
  # warning option is used if another warning is emitted in the same
  # compilation unit.
  # Therefore we invert the test for -Wno-... flags:
  #
  string(STRIP "${_flag}" _flag_stripped)
  set(_flag_sanitized "${_flag_stripped}")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    string(REPLACE "-Wno-" "-W" _flag_sanitized "${_flag_stripped}")
  endif()

  if(NOT "${_flag_stripped}" STREQUAL "")
    string(REGEX REPLACE "^-" "" _flag_name "${_flag_stripped}")
    string(REGEX REPLACE "\[-+,= \]" "_" _flag_name "${_flag_name}")

    CHECK_CXX_COMPILER_FLAG("${_flag_sanitized}" DEAL_II_HAVE_FLAG_${_flag_name})

    if(DEAL_II_HAVE_FLAG_${_flag_name})
      set(${_variable} "${${_variable}} ${_flag_stripped}")
      string(STRIP "${${_variable}}" ${_variable})
    endif()
  endif()

  reset_cmake_required()
endmacro()
