#
# A small macro used for (string-)appending a string flags to a
# string variable
#
# Usage:
#     ADD_FLAGS(variable flags)
#

MACRO(ADD_FLAGS variable flags)
  SET(${variable} "${${variable}} ${flags}")
ENDMACRO()

