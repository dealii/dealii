## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
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
# Replace all occurrences of "${flag}" with "${replacement}" in the string
# variable.
#
# Usage:
#     replace_flag(variable flag replacement)
#

macro(replace_flag _variable _flag _replacement)
  string(STRIP "${_replacement}" _replacement_stripped)
  string(REPLACE " " "  " ${_variable} "${${_variable}}")
  set(${_variable} " ${${_variable}} ")
  string(REPLACE " " "  " _flag2 "${_flag}")
  if(NOT "${_replacement_stripped}" STREQUAL "")
    string(REPLACE " ${_flag2} " " ${_replacement_stripped} " ${_variable} "${${_variable}}")
  else()
    string(REPLACE " ${_flag2} " " " ${_variable} "${${_variable}}")
  endif()
  string(REPLACE "  " " " ${_variable} "${${_variable}}")
  string(STRIP "${${_variable}}" ${_variable})
endmacro()
