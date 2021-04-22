## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2015 by the deal.II authors
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
# Search for a system library. In contrast to normal libraries we do this
# purely via "-l<library name>" instead of selecting a full library path..
#
# USAGE:
#   FIND_SYSTEM_LIBRARY(variable NAMES [list of possible names])
#

MACRO(FIND_SYSTEM_LIBRARY)
  SET(_argn ${ARGN})
  LIST(GET _argn 0 _variable)
  LIST(REMOVE_AT _argn 0 1)

  IF(NOT DEFINED ${_variable})
    FOREACH(_arg ${_argn})
      LIST(APPEND CMAKE_REQUIRED_LIBRARIES "${_arg}")
      CHECK_CXX_SOURCE_COMPILES("int main(){}" ${_variable})
      RESET_CMAKE_REQUIRED()

      IF(${_variable})
        UNSET(${_variable} CACHE)
        SET(${_variable} ${_arg} CACHE STRING "A system library.")
        SET(${_variable} ${_arg})
        MARK_AS_ADVANCED(${_variable})
        BREAK()
      ELSE()
        UNSET(${_variable} CACHE)
      ENDIF()
    ENDFOREACH()

    IF(NOT ${_variable})
      SET(${_variable} "${_variable}-NOTFOUND")
    ENDIF()
  ENDIF()
ENDMACRO()
