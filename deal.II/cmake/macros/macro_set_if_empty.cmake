#
# If 'variable' is empty it will be set to 'value'
#
MACRO(SET_IF_EMPTY variable value)
  IF("${${variable}}" STREQUAL "")
    SET(${variable} ${value})
  ENDIF()
ENDMACRO()

