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
# Remove all occurrences of "${flag}" in the string variable.
#
# Usage:
#     strip_flag(variable flag)
#

macro(strip_flag _variable _flag)
  string(REPLACE " " "  " ${_variable} "${${_variable}}")
  set(${_variable} " ${${_variable}} ")
  string(REPLACE " " "  " _flag2 "${_flag}")
  string(REPLACE " ${_flag2} " " " ${_variable} "${${_variable}}")
  string(REPLACE "  " " " ${_variable} "${${_variable}}")
  string(STRIP "${${_variable}}" ${_variable})
endmacro()

