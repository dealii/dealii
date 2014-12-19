## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
# A small macro used for converting a cmake list into a space
# separated string. This macro adds the string "prefix" in front of each
# element of the list.
#
# Usage:
#     TO_STRING_AND_ADD_PREFIX(string "prefix" ${list1} ${list2} ...)
#

MACRO(TO_STRING_AND_ADD_PREFIX _variable _prefix)
  SET(${_variable} "")
  FOREACH(_var ${ARGN})
    SET(${_variable} "${${_variable}} ${_prefix}${_var}")
  ENDFOREACH()
  STRING(STRIP "${${_variable}}" ${_variable})
ENDMACRO()
