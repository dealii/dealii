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
# Replace all occurences of "${flag}" with "${replacement}" in the string
# variable.
#
# Usage:
#     REPLACE_FLAG(variable flag replacement)
#

MACRO(REPLACE_FLAG _variable _flag _replacement)
  STRING(STRIP "${_replacement}" _replacement_stripped)
  STRING(REPLACE " " "  " ${_variable} "${${_variable}}")
  SET(${_variable} " ${${_variable}} ")
  STRING(REPLACE " " "  " _flag2 "${_flag}")
  IF(NOT "${_replacement_stripped}" STREQUAL "")
    STRING(REPLACE " ${_flag2} " " ${_replacement_stripped} " ${_variable} "${${_variable}}")
  ELSE()
    STRING(REPLACE " ${_flag2} " " " ${_variable} "${${_variable}}")
  ENDIF()
  STRING(REPLACE "  " " " ${_variable} "${${_variable}}")
  STRING(STRIP "${${_variable}}" ${_variable})
ENDMACRO()
