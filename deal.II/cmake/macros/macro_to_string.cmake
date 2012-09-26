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
