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
# A small macro to test whether a given list contains an element.
#
# Usage:
#     LIST_CONTAINS(var value list)
#
# var is set to true if list contains value as an element compared via
# STREQUAL.
#

MACRO(LIST_CONTAINS _var _value)
  SET(${_var})
  FOREACH (_value2 ${ARGN})
    IF("${_value}" STREQUAL "${_value2}")
      SET(${_var} TRUE)
      BREAK()
    ENDIF()
  ENDFOREACH()
ENDMACRO()

