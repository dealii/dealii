## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2021 by the deal.II authors
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
# A small wrapper around FIND_LIBRARY to be a bit more verbose
#

macro(deal_ii_find_library _library_name)
  # Save a string representation of the arguments before cmake's
  # FIND_FILE gets its hands on it.
  to_string(_str ${ARGN})

  find_library(${_library_name} ${ARGN})

  if(${_library_name} MATCHES "-NOTFOUND")
    message(STATUS "${_library_name} not found! The call was:")
    message(STATUS "    find_library(${_library_name} ${_str})")
  else()
    message(STATUS "Found ${_library_name}")
  endif()
endmacro()
