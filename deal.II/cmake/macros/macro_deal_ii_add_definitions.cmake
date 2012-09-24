
#
# TODO: A comment
#

MACRO(DEAL_II_ADD_DEFINITIONS name)

  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    GET_TARGET_PROPERTY(macro_definitions ${name}.g COMPILE_DEFINITIONS)
    SET_TARGET_PROPERTIES(${name}.g PROPERTIES
      COMPILE_DEFINITIONS "${ARGN};${macro_definitions}"
      )
  ENDIF()

  IF(CMAKE_BUILD_TYPE MATCHES "Release")
    GET_TARGET_PROPERTY(macro_definitions ${name} COMPILE_DEFINITIONS)
    SET_TARGET_PROPERTIES(${name} PROPERTIES
      COMPILE_DEFINITIONS "${ARGN};${macro_definitions}"
      )
  ENDIF()

ENDMACRO()

