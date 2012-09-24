
#
# TODO: A comment
#

MACRO(DEAL_II_ADD_C_LIBRARY library)

  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    #
    # and a debug target
    #
    ADD_LIBRARY(${library}_debug
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${library}_debug PROPERTIES
      LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_DEBUG}"
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_DEBUG}"
      COMPILE_FLAGS "${DEAL_II_C_FLAGS_DEBUG}"
      )

    FILE(APPEND
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/deal_ii_objects_debug
      "$<TARGET_OBJECTS:${library}_debug>\n"
      )
  ENDIF()

  IF(CMAKE_BUILD_TYPE MATCHES "Release")
    #
    # Add a release target
    #
    ADD_LIBRARY(${library}
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${library} PROPERTIES
      LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_RELEASE}"
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_RELEASE}"
      COMPILE_FLAGS "${DEAL_II_C_FLAGS_RELEASE}"
      )

    FILE(APPEND
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/deal_ii_objects
      "$<TARGET_OBJECTS:${library}>\n"
      )
  ENDIF()

ENDMACRO()
