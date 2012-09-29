
#
# Finalize the configuration:
#


#
# Hide some cmake specific cached variables. This is annoying...
#
MARK_AS_ADVANCED(file_cmd)


#
# Append the saved initial (cached) variables ${flags}_SAVED at the end of
# ${flags}, see setup_cached_compiler_flags.cmake and the main
# CMakeLists.txt for details.
#
FOREACH(flags ${deal_ii_used_flags})
  SET(${flags} "${${flags}} ${${flags}_SAVED}")
  #
  # Strip leading and trailing whitespace:
  #
  # STRING(STRIP "${flags}" flags)
  STRING(STRIP "${${flags}}" ${flags})
ENDFOREACH()


#
# Depulicate entries in DEAL_II_EXTERNAL_LIBRARIES(_...):
#
# Well, cmake has a bad habit of removing duplicate link library entries
# from targets. It wouldn't be that harmful if cmake would do this from
# behind... So, do it by hand to link correctly *sigh*
#
LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES)
LIST(REMOVE_DUPLICATES DEAL_II_EXTERNAL_LIBRARIES)
LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES)
FOREACH(build ${DEAL_II_BUILD_TYPES})
  LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES_${build})
  LIST(REMOVE_DUPLICATES DEAL_II_EXTERNAL_LIBRARIES_${build})
  LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES_${build})
ENDFOREACH()


#
# And print out a nice configuration summary:
#


MESSAGE("

*     *                                    *     *
*     *       deal.II configuration:       *     *
*     *                                    *     *


      CMAKE_BUILD_TYPE:       ${CMAKE_BUILD_TYPE}
      CMAKE_INSTALL_PREFIX:   ${CMAKE_INSTALL_PREFIX}
      CMAKE_SOURCE_DIR:       ${CMAKE_SOURCE_DIR}
      CMAKE_BINARY_DIR:       ${CMAKE_BINARY_DIR}

      CMAKE_CXX_COMPILER:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME}
                              ${CMAKE_CXX_COMPILER}

Compiler flags used for this build:
      CMAKE_CXX_FLAGS:                     ${CMAKE_CXX_FLAGS}")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  MESSAGE("      DEAL_II_CXX_FLAGS_RELEASE:           ${DEAL_II_CXX_FLAGS_RELEASE}")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  MESSAGE("      DEAL_II_CXX_FLAGS_DEBUG:             ${DEAL_II_CXX_FLAGS_DEBUG}")
ENDIF()
MESSAGE("      CMAKE_SHARED_LINKER_FLAGS:           ${CMAKE_SHARED_LINKER_FLAGS}")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  MESSAGE("      DEAL_II_SHARED_LINKER_FLAGS_RELEASE:  ${DEAL_II_SHARED_LINKER_FLAGS_RELEASE}")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  MESSAGE("      DEAL_II_SHARED_LINKER_FLAGS_DEBUG:   ${DEAL_II_SHARED_LINKER_FLAGS_DEBUG}")
ENDIF()


IF(FEATURE_UMFPACK_CONTRIB_CONFIGURED)
  MESSAGE("
The contrib UMFPACK library will be compiled with the following C compiler:
      CMAKE_C_COMPILER:         ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}
                                ${CMAKE_C_COMPILER}
      CMAKE_C_FLAGS:           ${CMAKE_C_FLAGS}")
  IF(CMAKE_BUILD_TYPE MATCHES "Release")
    MESSAGE("      DEAL_II_C_FLAGS_RELEASE: ${DEAL_II_C_FLAGS_RELEASE}")
  ENDIF()
  IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    MESSAGE("      DEAL_II_C_FLAGS_DEBUG:   ${DEAL_II_C_FLAGS_DEBUG}")
  ENDIF()
ENDIF()

IF(NOT DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
  MESSAGE("\n"
    "WARNING: DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS is set to OFF\n"
    )
ELSE()
  IF(NOT DEAL_II_KNOWN_COMPILER)
    MESSAGE("\n"
      "WARNING: Unknown compiler! Please set compiler flags by hand.\n"
      )
  ENDIF()
ENDIF()

MESSAGE("
  Configured Features ("
  "DEAL_II_FEATURE_AUTODETECTION = ${DEAL_II_FEATURE_AUTODETECTION}, "
  "DEAL_II_ALLOW_CONTRIB = ${DEAL_II_ALLOW_CONTRIB}):"
  )
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
        IF(DEAL_II_FORCE_CONTRIB_${feature})
          MESSAGE("      ${var} set up with contrib packages (forced)")
        ELSE()
          MESSAGE("      ${var} set up with contrib packages")
        ENDIF()
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
      MESSAGE("      ${var}")
    ELSE()
      MESSAGE("    ( ${var} = ${${var}} )")
    ENDIF()
  ENDIF()
ENDFOREACH()


MESSAGE("\n")

