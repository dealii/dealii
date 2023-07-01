## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2022 by the deal.II authors
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

macro(deal_ii_find_program _file_name)
  # Save a string representation of the arguments before cmake's
  # FIND_PROGRAM gets its hands on it.
  to_string(_str ${ARGN})

  find_program(${_file_name} ${ARGN})

  if(${_file_name} MATCHES "-NOTFOUND")
    message(STATUS "${_file_name} not found! The call was:")
    message(STATUS "    find_program(${_file_name} ${_str})")
  else()
    message(STATUS "Found ${_file_name}")
  endif()
endmacro()
