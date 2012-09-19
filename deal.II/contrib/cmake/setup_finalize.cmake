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

   *                                    *
   *  deal.II successfully configured!  *
   *                                    *


      CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}
      CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}
      CMAKE_SOURCE_DIR:     ${CMAKE_SOURCE_DIR}
      CMAKE_BINARY_DIR:     ${CMAKE_BINARY_DIR}

      CMAKE_CXX_COMPILER:   ${CMAKE_CXX_COMPILER_ID} "
"${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME}

General compiler flags (used by all build targets):

      CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}
")

IF(CMAKE_BUILD_TYPE MATCHES "Release")
  MESSAGE("Additional compiler flags used for the Release target:

      CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}
")
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  MESSAGE("Additional compiler flags used for the Debug target:

      CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}
")
ENDIF()

MESSAGE("Configured linker flags:

      CMAKE_SHARED_LINKER_FLAGS: ${CMAKE_SHARED_LINKER_FLAGS}
")

IF(FEATURE_UMFPACK_CONTRIB_CONFIGURED)
  MESSAGE("The contrib UMFPACK library will be compiled with the following C compiler:

      CMAKE_C_COMPILER:       ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}
      CMAKE_C_FLAGS:         ${CMAKE_C_FLAGS}
      CMAKE_C_FLAGS_DEBUG:   ${CMAKE_C_FLAGS_DEBUG}
      CMAKE_C_FLAGS_RELEASE: ${CMAKE_C_FLAGS_RELEASE}

")
ENDIF()

MESSAGE("
Configured Features (DEAL_II_FEATURE_AUTODETECTION = ${DEAL_II_FEATURE_AUTODETECTION}):
")

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

")
