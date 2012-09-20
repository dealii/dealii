#
# A small macro to test whether a given list contains an element.
#
# Usage:
#     LIST_CONTAINS(var value list)
#
# var is set to true if list contains value as an element compared via
# STREQUAL.
#

MACRO(LIST_CONTAINS var value)
  SET(${var})
  FOREACH (value2 ${ARGN})
    IF (${value} STREQUAL ${value2})
      SET(${var} TRUE)
    ENDIF (${value} STREQUAL ${value2})
  ENDFOREACH (value2)
ENDMACRO(LIST_CONTAINS)

