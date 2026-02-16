## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2012 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# A small macro to test whether a given list contains an element.
#
# Usage:
#     item_matches(var regex list)
#
# var is set to true if list contains an item that matches regex.
#

macro(item_matches _var _regex)
  set(${_var})
  foreach (_item ${ARGN})
    if("${_item}" MATCHES ${_regex})
      set(${_var} TRUE)
      break()
    endif()
  endforeach()
endmacro()
