## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2024 by the deal.II authors
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
