
#
# A small macro used for converting a cmake list <list> into a space
# separated <string>. This macro adds the string <prefix> in front of each
# element of the list <list>.
#
# Usage:
#     TO_STRING_AND_ADD_PREFIX(string list)
#

MACRO(TO_STRING_AND_ADD_PREFIX variable prefix list)
  SET(${variable} "")
  FOREACH(var ${${list}})
    SET(${variable} "${${variable}} ${prefix}${var}")
  ENDFOREACH()
  STRING(STRIP "${${variable}}" ${variable})
ENDMACRO()
