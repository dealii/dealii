#
# Remove the first occurence of flag in the string variable.
#
# Usage:
#     STRIP_FLAG(variable flag)
#

MACRO(STRIP_FLAG variable flag)
  IF(NOT "${variable}" STREQUAL "")
    STRING(REGEX REPLACE " ${flag}" "" ${variable} ${${variable}})
  ENDIF()
ENDMACRO()

