
#
# TODO: A comment
#

MACRO(DEAL_II_ADD_DEPENDENCIES name target)

  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    ADD_DEPENDENCIES(${name}.g ${target}.g)
  ENDIF()

  IF(CMAKE_BUILD_TYPE MATCHES "Release")
    ADD_DEPENDENCIES(${name} ${target})
  ENDIF()

ENDMACRO()
