
MACRO(SET_IF_EMPTY variable value)
  IF(NOT ${${variable}} STREQUAL "")
    SET(${variable} ${value})
ENDMACRO()

