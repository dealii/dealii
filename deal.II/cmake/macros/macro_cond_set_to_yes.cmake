
#
# IF(bool), set variable to "yes".
#
# Usage:
#     COND_SET_TO_YES(bool variable)
#

MACRO(COND_SET_TO_YES bool variable)
  IF(${bool})
    SET(${variable} "yes")
  ELSE()
    SET(${variable} "no")
  ENDIF()
ENDMACRO()

