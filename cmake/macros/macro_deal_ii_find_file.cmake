## ---------------------------------------------------------------------
##
## Copyright (C) 2014 by the deal.II authors
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
# A small wrapper around FIND_FILE to be a bit more verbose
#

MACRO(DEAL_II_FIND_FILE _file_name)
  FIND_FILE(${_file_name} ${ARGN})

  IF(${_file_name} MATCHES "-NOTFOUND")
    MESSAGE(STATUS "${_file_name} not found! Call:")
    TO_STRING(_str ${ARGN})
    MESSAGE(STATUS "    FIND_FILE(${_file_name} ${_str})")
  ELSE()
    MESSAGE(STATUS "Found ${_file_name}")
  ENDIF()
ENDMACRO()
