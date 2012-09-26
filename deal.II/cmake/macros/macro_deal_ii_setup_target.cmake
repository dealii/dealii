#
# This file implements the DEAL_II_SETUP_TARGET macro, which is
# part of the deal.II library.
#
# Usage:
#       DEAL_II_SETUP_TARGET(target)
#
# This macro sets the necessary include directories, linker flags, compile
# definitions and the external libraries the target will be linked against.
#

MACRO(DEAL_II_SETUP_TARGET target)

  CMAKE_MINIMUM_REQUIRED(2.8.8)

  SET_TARGET_PROPERTIES(${target} PROPERTIES
    INCLUDE_DIRECTORIES
      "${DEAL_II_EXTERNAL_INCLUDE_DIRS};${DEAL_II_INCLUDE_DIRS}"
    LINK_FLAGS
      "${DEAL_II_LINKER_FLAGS}"
    LINK_FLAGS_DEBUG
      "${DEAL_II_LINKER_FLAGS_DEBUG}"
    LINK_FLAGS_RELEASE
      "${DEAL_II_LINKER_FLAGS_RELEASE}"
    COMPILE_DEFINITIONS
      "${DEAL_II_USER_DEFINITIONS}"
    COMPILE_DEFINITIONS_DEBUG
      "${DEAL_II_USER_DEFINITIONS_DEBUG}"
    COMPILE_DEFINITIONS_RELEASE
      "${DEAL_II_USER_DEFINITIONS_RELEASE}"
    )

  TARGET_LINK_LIBRARIES(${target}
    debug
      ${DEAL_II_LIBRARIES_DEBUG}
      ${DEAL_II_EXTERNAL_LIBRARIES_DEBUG}
    optimized
      ${DEAL_II_LIBRARIES_RELEASE}
      ${DEAL_II_EXTERNAL_LIBRARIES_RELEASE}
    general
      ${DEAL_II_EXTERNAL_LIBRARIES}
    )

ENDMACRO()

