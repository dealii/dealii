#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

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
IF(NOT "${DEAL_II_EXTERNAL_LIBRARIES}" STREQUAL "")
  LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES)
  LIST(REMOVE_DUPLICATES DEAL_II_EXTERNAL_LIBRARIES)
  LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES)
ENDIF()
FOREACH(build ${DEAL_II_BUILD_TYPES})
  IF(NOT "${DEAL_II_EXTERNAL_LIBRARIES_${build}}" STREQUAL "")
    LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES_${build})
    LIST(REMOVE_DUPLICATES DEAL_II_EXTERNAL_LIBRARIES_${build})
    LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES_${build})
  ENDIF()
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


IF(FEATURE_UMFPACK_BUNDLED_CONFIGURED)
  MESSAGE("
The bundled UMFPACK library will be compiled with the following C compiler:
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

MESSAGE("\nConfigured Features ("
  "DEAL_II_FEATURE_AUTODETECTION = ${DEAL_II_FEATURE_AUTODETECTION}, "
  "DEAL_II_ALLOW_BUNDLED = ${DEAL_II_ALLOW_BUNDLED}):"
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

      IF(FEATURE_${feature}_BUNDLED_CONFIGURED)
        IF(DEAL_II_FORCE_BUNDLED_${feature})
          MESSAGE("      ${var} set up with bundled packages (forced)")
        ELSE()
          MESSAGE("      ${var} set up with bundled packages")
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

