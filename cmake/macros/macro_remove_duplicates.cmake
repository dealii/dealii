## ---------------------------------------------------------------------
##
## Copyright (C) 2014 by the deal.II authors
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
