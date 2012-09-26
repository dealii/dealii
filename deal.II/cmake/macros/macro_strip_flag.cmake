#
# Remove all occurences of "${flag}" in the string variable.
#
# Usage:
#     STRIP_FLAG(variable flag)
#

MACRO(STRIP_FLAG variable flag)
  SET(${variable} " ${${variable}}")
  STRING(REPLACE " ${flag}" "" "${variable}" ${${variable}})
  STRING(STRIP "${${variable}}" ${variable})
ENDMACRO()

