
#
# TODO: A comment
#

MACRO(DEAL_II_ADD_C_LIBRARY library)

  FOREACH(build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${build} build_lowercase)

    ADD_LIBRARY(${library}.${build_lowercase}
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${library}.${build_lowercase} PROPERTIES
      LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_${build}}"
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${build}}"
      COMPILE_FLAGS "${DEAL_II_C_FLAGS_${build}}"
      )

    FILE(APPEND
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_objects_${build_lowercase}
      "$<TARGET_OBJECTS:${library}.${build_lowercase}>\n"
      )
  ENDFOREACH()

ENDMACRO()
