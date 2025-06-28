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
# Remove all occurrences of "${flag}" in the string variable.
#
# Usage:
#     strip_flag(variable flag)
#

macro(strip_flag _variable _flag)
  string(REPLACE " " "  " ${_variable} "${${_variable}}")
  set(${_variable} " ${${_variable}} ")
  string(REPLACE " " "  " _flag2 "${_flag}")
  string(REPLACE " ${_flag2} " " " ${_variable} "${${_variable}}")
  string(REPLACE "  " " " ${_variable} "${${_variable}}")
  string(STRIP "${${_variable}}" ${_variable})
endmacro()
