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
# Check for a compiler bug.
#
# Usage:
#     CHECK_CXX_COMPILER_BUG(source var),
#
# where source is a snipped of source code and var is a variable that will
# be set to true if the source could not be compiled and linked successfully.
# (This just inverts the logic of CHECK_CXX_SOURCE_COMPILES.)
#

MACRO(CHECK_CXX_COMPILER_BUG _source _var)
  IF(NOT DEFINED ${_var}_OK)
    CHECK_CXX_SOURCE_COMPILES(
      "${_source}"
      ${_var}_OK
      )
    IF(${_var}_OK)
      MESSAGE(STATUS "Test successful, do not define ${_var}")
    ELSE()
      MESSAGE(STATUS "Test unsuccessful, define ${_var}")
    ENDIF()
  ENDIF()

  IF(${_var}_OK)
    SET(${_var})
  ELSE()
    SET(${_var} TRUE)
  ENDIF()
ENDMACRO()

