
#
# A small macro used for converting a cmake list <list> into a space
# separated <string>:
#
# Usage:
#     TO_STRING(string list)
#

MACRO(TO_STRING variable list)
  SET(${variable} "")
  FOREACH(var ${${list}})
    SET(${variable} "${${variable}} ${var}")
  ENDFOREACH()
  STRING(STRIP "${${variable}}" ${variable})
ENDMACRO()
