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
# Tests whether it is possible to compile and link a dummy program with a
# given flag. If so, add it to variable.
#
# Note: This macro will reset the CMAKE_REQUIRED_* variables.
#
# Usage:
#     enable_if_links(variable flag)
#

macro(enable_if_links _variable _flag)
  reset_cmake_required()
  enable_if_supported(CMAKE_REQUIRED_FLAGS "-Werror")

  string(STRIP "${_flag}" _flag_stripped)
  if(NOT "${_flag_stripped}" STREQUAL "")
    string(REGEX REPLACE "^-" "" _flag_name "${_flag_stripped}")
    string(REGEX REPLACE "\[-+,\]" "_" _flag_name "${_flag_name}")

    list(APPEND CMAKE_REQUIRED_LIBRARIES "${_flag_stripped}")
    CHECK_CXX_SOURCE_COMPILES("int main(){}" DEAL_II_HAVE_LINKER_FLAG_${_flag_name})

    if(DEAL_II_HAVE_LINKER_FLAG_${_flag_name})
      set(${_variable} "${${_variable}} ${_flag_stripped}")
      string(STRIP "${${_variable}}" ${_variable})
    endif()
  endif()

  reset_cmake_required()
endmacro()
