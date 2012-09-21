
#
# A small macro used for converting a cmake list into a space
# separated string. This macro adds the string "prefix" in front of each
# element of the list.
#
# Usage:
#     TO_STRING_AND_ADD_PREFIX(string "prefix" ${list1} ${list2} ...)
#

MACRO(TO_STRING_AND_ADD_PREFIX variable prefix)
  SET(${variable} "")
  FOREACH(var ${ARGN})
    SET(${variable} "${${variable}} ${prefix}${var}")
  ENDFOREACH()
  STRING(STRIP "${${variable}}" ${variable})
ENDMACRO()
