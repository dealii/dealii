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
# A small macro used for converting a list into a space
# separated string:
#
# Usage:
#     TO_STRING(string ${list1} ${list2} ...)
#

MACRO(TO_STRING variable)
  SET(${variable} "")
  FOREACH(var  ${ARGN})
    SET(${variable} "${${variable}} ${var}")
  ENDFOREACH()
  STRING(STRIP "${${variable}}" ${variable})
ENDMACRO()
