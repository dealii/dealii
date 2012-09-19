
MACRO(STRIP_FLAG variable flag)
  IF(NOT "${variable}" STREQUAL "")
    STRING(REGEX REPLACE " ${flag}" "" ${variable} ${${variable}})
  ENDIF()
ENDMACRO()

