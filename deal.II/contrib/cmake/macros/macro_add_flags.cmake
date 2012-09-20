#
# A small macro used for (string-)appending flags to a variable
#
# Usage:
#     ADD_FLAGS(variable flag)
#

MACRO(ADD_FLAGS variable flag)
  SET(${variable} "${${variable}} ${flag}")
ENDMACRO()

