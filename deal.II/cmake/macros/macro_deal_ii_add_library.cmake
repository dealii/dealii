#
# A small wrapper around ADD_LIBRARY that will define a target for each
# build type specified in DEAL_II_BUILD_TYPES
#
# It is assumed that the desired compilation configuration is set via
#   DEAL_II_SHARED_LINKER_FLAGS_${build}
#   DEAL_II_CXX_FLAGS_${build}
#   DEAL_II_DEFINITIONS_${build}
#
# as well as the global (for all build types)
#   CMAKE_SHARED_LINKER_FLAGS
#   CMAKE_CXX_FLAGS
#   DEAL_II_DEFINITIONS
#

MACRO(DEAL_II_ADD_LIBRARY library)

  FOREACH(build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${build} build_lowercase)

    ADD_LIBRARY(${library}.${build_lowercase}
      ${ARGN}
      )

    SET_TARGET_PROPERTIES(${library}.${build_lowercase} PROPERTIES
      LINK_FLAGS "${DEAL_II_SHARED_LINKER_FLAGS_${build}}"
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${build}}"
      COMPILE_FLAGS "${DEAL_II_CXX_FLAGS_${build}}"
      )

    FILE(APPEND
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_objects_${build_lowercase}
      "$<TARGET_OBJECTS:${library}.${build_lowercase}>\n"
      )
  ENDFOREACH()

ENDMACRO()
