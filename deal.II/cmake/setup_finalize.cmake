#
# Finalize the configuration:
#


#
# Dependency check:
# DEAL_II_COMPONENT_DOCUMENTATION needs DEAL_II_WITH_DOXYGEN.
# TODO: It is a bit sloppy to test this here. But this is the only
# dependency of this kind atm.
#
IF(DEAL_II_COMPONENT_DOCUMENTATION AND NOT DEAL_II_WITH_DOXYGEN)
  MESSAGE(SEND_ERROR "\n"
    "DEAL_II_COMPONENT_DOCUMENTATION has unmet configuration requirements: "
    "DEAL_II_WITH_DOXYGEN has to be set to \"ON\".\n\n"
    )
ENDIF()



#
# Append the saved initial (cached) variables ${flags}_SAVED at the end of
# ${flags}, see setup_cached_compiler_flags.cmake and the main
# CMakeLists.txt for details.
#
FOREACH(flags ${deal_ii_used_flags})
  SET(${flags} "${${flags}} ${${flags}_SAVED}")
ENDFOREACH()



#
# And print out a nice configuration summary:
#


MESSAGE("

*     *                                    *     *
*     *  deal.II successfully configured!  *     *
*     *                                    *     *


      CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}
      CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}
      CMAKE_SOURCE_DIR:     ${CMAKE_SOURCE_DIR}
      CMAKE_BINARY_DIR:     ${CMAKE_BINARY_DIR}

      CMAKE_CXX_COMPILER:   ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME}

Compiler flags used for this build:
      CMAKE_CXX_FLAGS:           ${CMAKE_CXX_FLAGS}")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  MESSAGE("      CMAKE_CXX_FLAGS_RELEASE:   ${CMAKE_CXX_FLAGS_RELEASE}")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  MESSAGE("      CMAKE_CXX_FLAGS_DEBUG:     ${CMAKE_CXX_FLAGS_DEBUG}")
ENDIF()
MESSAGE("      CMAKE_SHARED_LINKER_FLAGS: ${CMAKE_SHARED_LINKER_FLAGS}")


IF(FEATURE_UMFPACK_CONTRIB_CONFIGURED)
  MESSAGE("
The contrib UMFPACK library will be compiled with the following C compiler:
      CMAKE_C_COMPILER:       ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}
      CMAKE_C_FLAGS:         ${CMAKE_C_FLAGS}")
  IF(CMAKE_BUILD_TYPE MATCHES "Release")
    MESSAGE("      CMAKE_C_FLAGS_RELEASE: ${CMAKE_C_FLAGS_RELEASE}")
  ENDIF()
  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    MESSAGE("      CMAKE_C_FLAGS_DEBUG:   ${CMAKE_C_FLAGS_DEBUG}")
  ENDIF()
ENDIF()


IF(NOT DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
  MESSAGE("\n
WARNING: DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS is set to OFF\n")
ENDIF()


MESSAGE("
Configured Features (DEAL_II_FEATURE_AUTODETECTION = ${DEAL_II_FEATURE_AUTODETECTION}):")
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  IF(var MATCHES "DEAL_II_WITH")
    IF(${${var}})
      # FEATURE is enabled
      STRING(REGEX REPLACE "^DEAL_II_WITH_" "" feature ${var})

      IF(FEATURE_${feature}_EXTERNAL_CONFIGURED)
        MESSAGE("      ${var} set up with external dependencies")
      ENDIF()

      IF(FEATURE_${feature}_CONTRIB_CONFIGURED)
        MESSAGE("      ${var} set up with contrib packages")
      ENDIF()
    ELSE()
      # FEATURE is disabled
      MESSAGE("    ( ${var} = ${${var}} )")
    ENDIF()
  ENDIF()
ENDFOREACH()


MESSAGE("
Component configuration:")
FOREACH(var ${res})
  IF(var MATCHES "DEAL_II_COMPONENT")
    IF(${${var}})
      MESSAGE("      ${var} = ${${var}}")
    ELSE()
      MESSAGE("    ( ${var} = ${${var}} )")
    ENDIF()
  ENDIF()
ENDFOREACH()

MESSAGE("
Further configuration:")
FOREACH(var ${res})
  IF(var MATCHES "DEAL_II_INSTALL")
    IF(${${var}})
      MESSAGE("      ${var} = ${${var}}")
    ELSE()
      MESSAGE("    ( ${var} = ${${var}} )")
    ENDIF()
  ENDIF()
ENDFOREACH()


MESSAGE("\n")
