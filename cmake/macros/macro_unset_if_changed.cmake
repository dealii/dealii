## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
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
# Usage:
#   UNSET_IF_CHANGED(<internal variable> "string"
#     [cached variable names]
#     )
#
# This macro caches the supplied "string" internally in ${<internal
# variable>} and unsets all supplied (cached) variables if this string
# changes.
#
MACRO(UNSET_IF_CHANGED _variable _string)
  IF(DEFINED ${_variable})
    IF(NOT "${${_variable}}" STREQUAL "${_string}")
      FOREACH(_arg ${ARGN})
        MESSAGE(STATUS
          "Configuration changed. Unsetting cached variable \"${_arg}\" and rerunning checks.")
        UNSET(${_arg} CACHE)
      ENDFOREACH()
    ENDIF()
  ENDIF()
  SET(${_variable} "${_string}" CACHE INTERNAL "" FORCE)
ENDMACRO()
