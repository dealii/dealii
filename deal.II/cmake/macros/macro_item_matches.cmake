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
#     ITEM_MATCHES(var regex list)
#
# var is set to true if list contains an item that matches regex.
#

MACRO(ITEM_MATCHES _var _regex)
  SET(${_var})
  FOREACH (_item ${ARGN})
    IF("${_item}" MATCHES ${_regex})
      SET(${_var} TRUE)
      BREAK()
    ENDIF()
  ENDFOREACH()
ENDMACRO()

