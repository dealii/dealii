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
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Remove duplicate entries from a list. Optionally do this in reverse
# order, keeping the rightmost element
#
# Usage:
#     REMOVE_DUPLICATES(list [REVERSE])
#

MACRO(REMOVE_DUPLICATES _list)
  IF(NOT "${${_list}}" STREQUAL "")
    IF("${ARGN}" STREQUAL "REVERSE")
      LIST(REVERSE ${_list})
    ENDIF()
    LIST(REMOVE_DUPLICATES ${_list})
    IF("${ARGN}" STREQUAL "REVERSE")
      LIST(REVERSE ${_list})
    ENDIF()
  ENDIF()
ENDMACRO()
