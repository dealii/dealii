
#
# TODO: A comment
#

MACRO(DEAL_II_ADD_C_LIBRARY name)

  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    #
    # and a debug target
    #
    ADD_LIBRARY(${name}.g
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${name}.g PROPERTIES
      LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_DEBUG}"
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_DEBUG}"
      COMPILE_FLAGS "${DEAL_II_C_FLAGS_DEBUG}"
      )

    SET(deal_ii_objects.g
      ${deal_ii_objects.g}
      $<TARGET_OBJECTS:${name}.g>
      PARENT_SCOPE
      )

  ENDIF()

  IF(CMAKE_BUILD_TYPE MATCHES "Release")
    #
    # Add a release target
    #
    ADD_LIBRARY(${name}
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${name} PROPERTIES
      LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_RELEASE}"
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_RELEASE}"
      COMPILE_FLAGS "${DEAL_II_C_FLAGS_RELEASE}"
      )

    SET(deal_ii_objects
      ${deal_ii_objects}
      $<TARGET_OBJECTS:${name}>
      PARENT_SCOPE
      )
  ENDIF()

ENDMACRO()
