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
# Replace all occurences of "${flag}" with "${replacement}" in the string
# variable.
#
# Usage:
#     STRIP_FLAG(variable flag replacement)
#

MACRO(REPLACE_FLAG _variable _flag _replacement)
  STRING(STRIP "${_replacement}" _replacement_stripped)
  STRING(REPLACE " " "  " ${_variable} "${${_variable}}")
  SET(${_variable} " ${${_variable}} ")
  STRING(REPLACE " " "  " _flag2 "${_flag}")
  IF(NOT "${_replacement_stripped}" STREQUAL "")
    STRING(REPLACE " ${_flag2} " " ${_replacement_stripped} " ${_variable} "${${_variable}}")
  ELSE()
    STRING(REPLACE " ${_flag2} " " " ${_variable} "${${_variable}}")
  ENDIF()
  STRING(REPLACE "  " " " ${_variable} "${${_variable}}")
  STRING(STRIP "${${_variable}}" ${_variable})
ENDMACRO()
