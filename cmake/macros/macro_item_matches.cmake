## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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

