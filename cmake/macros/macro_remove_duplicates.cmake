## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2022 by the deal.II authors
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
# Remove duplicate entries from a list. Optionally do this in reverse
# order, keeping the rightmost element
#
# Usage:
#     remove_duplicates(list [REVERSE])
#

macro(remove_duplicates _list)
  if(NOT "${${_list}}" STREQUAL "")
    if("${ARGN}" STREQUAL "REVERSE")
      list(REVERSE ${_list})
    endif()
    list(REMOVE_DUPLICATES ${_list})
    if("${ARGN}" STREQUAL "REVERSE")
      list(REVERSE ${_list})
    endif()
  endif()
endmacro()
