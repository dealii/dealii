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
# A small macro used for (string-)appending a string "${flags}" to a
# string "${variable}"
#
# Usage:
#     add_flags(variable flags)
#

macro(add_flags _variable _flags)
  string(STRIP "${_flags}" _flags_stripped)
  if(NOT "${_flags_stripped}" STREQUAL "")
    set(${_variable} "${${_variable}} ${_flags}")
    string(STRIP "${${_variable}}" ${_variable})
  endif()
endmacro()

