
#
# TODO: A comment
#

MACRO(DEAL_II_ADD_DEPENDENCIES name target)

  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    ADD_DEPENDENCIES(${name}_debug ${target}_debug)
  ENDIF()

  IF(CMAKE_BUILD_TYPE MATCHES "Release")
    ADD_DEPENDENCIES(${name} ${target})
  ENDIF()

ENDMACRO()
