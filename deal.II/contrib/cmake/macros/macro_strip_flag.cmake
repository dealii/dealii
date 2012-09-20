#
# Remove all occurences of "${flag}" in the string variable.
#
# Usage:
#     STRIP_FLAG(variable flag)
#

MACRO(STRIP_FLAG variable flag)
  IF(NOT "${variable}" STREQUAL "")
    SET(${variable} " ${${variable}}")
    STRING(REPLACE " ${flag}" "" ${variable} ${${variable}})
    IF(NOT "${variable}" STREQUAL "")
      STRING(STRIP ${${variable}} ${variable})
    ENDIF()
  ENDIF()
ENDMACRO()

