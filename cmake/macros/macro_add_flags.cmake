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
# A small macro used for (string-)appending a string "${flags}" to a
# string "${variable}"
#
# Usage:
#     ADD_FLAGS(variable flags)
#

MACRO(ADD_FLAGS _variable _flags)
  STRING(STRIP "${_flags}" _flags_stripped)
  IF(NOT "${_flags_stripped}" STREQUAL "")
    SET(${_variable} "${${_variable}} ${_flags}")
    STRING(STRIP "${${_variable}}" ${_variable})
  ENDIF()
ENDMACRO()

