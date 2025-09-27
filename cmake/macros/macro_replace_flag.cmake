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
# Replace all occurrences of "${flag}" with "${replacement}" in the string
# variable.
#
# Usage:
#     replace_flag(variable flag replacement)
#

macro(replace_flag _variable _flag _replacement)
  string(STRIP "${_replacement}" _replacement_stripped)
  string(REPLACE " " "  " ${_variable} "${${_variable}}")
  set(${_variable} " ${${_variable}} ")
  string(REPLACE " " "  " _flag2 "${_flag}")
  if(NOT "${_replacement_stripped}" STREQUAL "")
    string(REPLACE " ${_flag2} " " ${_replacement_stripped} " ${_variable} "${${_variable}}")
  else()
    string(REPLACE " ${_flag2} " " " ${_variable} "${${_variable}}")
  endif()
  string(REPLACE "  " " " ${_variable} "${${_variable}}")
  string(STRIP "${${_variable}}" ${_variable})
endmacro()
