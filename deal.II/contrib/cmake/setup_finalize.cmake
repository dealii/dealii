#
# Append the saved cache variable ${flags}_SAVED at the end of ${flags}
#
FOREACH(flags ${deal_ii_used_flags})
  SET(${flags} "${${flags}} ${${flags}_SAVED}")
ENDFOREACH()

#
# And print out a nice configuration summary:
#
MESSAGE("

*                                        *
*    deal.II successfully configured!    *
*                                        *


    CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}
    CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}

General compiler flags (used by all build targets):

    CMAKE_C_FLAGS:       ${CMAKE_C_FLAGS}
    CMAKE_CXX_FLAGS:     ${CMAKE_CXX_FLAGS}
")

IF(CMAKE_BUILD_TYPE MATCHES "Release")
  MESSAGE("
Additional compiler flags used for the Release target:

    CMAKE_C_FLAGS_RELEASE:   ${CMAKE_C_FLAGS_RELEASE}
    CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}
")
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  MESSAGE("
Additional compiler flags used for the Debug target:

    CMAKE_C_FLAGS_DEBUG:   ${CMAKE_C_FLAGS_DEBUG}
    CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}
")
ENDIF()

MESSAGE("
(Note: Flags set with ccmake or the command line will be appended at the end of the default configuration)


Configured Features:
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
