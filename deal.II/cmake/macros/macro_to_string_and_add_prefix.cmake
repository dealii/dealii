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
# A small macro used for converting a cmake list into a space
# separated string. This macro adds the string "prefix" in front of each
# element of the list.
#
# Usage:
#     TO_STRING_AND_ADD_PREFIX(string "prefix" ${list1} ${list2} ...)
#

MACRO(TO_STRING_AND_ADD_PREFIX _variable _prefix)
  SET(${_variable} "")
  FOREACH(_var ${ARGN})
    SET(${_variable} "${${_variable}} ${_prefix}${_var}")
  ENDFOREACH()
  STRING(STRIP "${${_variable}}" ${_variable})
ENDMACRO()
