## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the deal.II authors
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
# A small wrapper around FIND_PROGRAM to be a bit more verbose
#

MACRO(DEAL_II_FIND_PROGRAM _file_name)
  # Save a string representation of the arguments before cmake's
  # FIND_PROGRAM gets its hands on it.
  TO_STRING(_str ${ARGN})

  FIND_PROGRAM(${_file_name} ${ARGN})

  IF(${_file_name} MATCHES "-NOTFOUND")
    MESSAGE(STATUS "${_file_name} not found! The call was:")
    MESSAGE(STATUS "    FIND_PROGRAM(${_file_name} ${_str})")
  ELSE()
    MESSAGE(STATUS "Found ${_file_name}")
  ENDIF()
ENDMACRO()
