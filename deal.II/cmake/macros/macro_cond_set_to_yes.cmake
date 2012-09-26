#
# If bool is "true" (in a cmake fashion...), set variable to "yes",
# otherwise to "no".
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

