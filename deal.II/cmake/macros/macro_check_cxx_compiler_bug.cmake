#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

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
  CHECK_CXX_SOURCE_COMPILES(
    "${_source}"
    ${_var}_OK)

  IF(${_var}_OK)
    MESSAGE(STATUS "Test successful, do not define ${_var}")
  ELSE()
    MESSAGE(STATUS "Test unsuccessful, define ${_var}")
    SET(${_var} TRUE)
  ENDIF()
ENDMACRO()

