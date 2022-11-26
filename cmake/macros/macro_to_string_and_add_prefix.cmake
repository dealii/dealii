## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
# A small macro used for converting a cmake list into a space
# separated string. This macro adds the string "prefix" in front of each
# element of the list.
#
# Usage:
#     to_string_and_add_prefix(string "prefix" ${list1} ${list2} ...)
#

macro(to_string_and_add_prefix _variable _prefix)
  set(${_variable} "")
  foreach(_var ${ARGN})
    set(${_variable} "${${_variable}} ${_prefix}${_var}")
  endforeach()
  string(STRIP "${${_variable}}" ${_variable})
endmacro()
