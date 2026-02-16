## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# A small macro used for converting a list into a space
# separated string:
#
# Usage:
#     to_string(string ${list1} ${list2} ...)
#

macro(to_string _variable)
  set(${_variable} "")
  foreach(_var  ${ARGN})
    set(${_variable} "${${_variable}} ${_var}")
  endforeach()
  string(STRIP "${${_variable}}" ${_variable})
endmacro()
