#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
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


###########################################################################
#                                                                         #
#                       Finalize the configuration:                       #
#                                                                         #
###########################################################################

#
# Hide some cmake specific cached variables. This is annoying...
#
MARK_AS_ADVANCED(file_cmd)

#
# Append the saved initial (cached) variables ${flags}_SAVED at the end of
# ${flags}, see setup_cached_compiler_flags.cmake and the main
# CMakeLists.txt for details.
#
FOREACH(_flags ${DEAL_II_USED_FLAGS})
  # Strip leading and trailing whitespace:
  STRING(STRIP "${${_flags}} ${${_flags}_SAVED}" ${_flags})
ENDFOREACH()

#
# Depulicate entries in DEAL_II_EXTERNAL_LIBRARIES(_...):
#
IF(NOT "${DEAL_II_EXTERNAL_LIBRARIES}" STREQUAL "")
  LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES)
  LIST(REMOVE_DUPLICATES DEAL_II_EXTERNAL_LIBRARIES)
  LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES)
ENDIF()
FOREACH(_build ${DEAL_II_BUILD_TYPES})
  IF(NOT "${DEAL_II_EXTERNAL_LIBRARIES_${_build}}" STREQUAL "")
    LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES_${_build})
    LIST(REMOVE_DUPLICATES DEAL_II_EXTERNAL_LIBRARIES_${_build})
    LIST(REVERSE DEAL_II_EXTERNAL_LIBRARIES_${_build})
  ENDIF()
ENDFOREACH()


###########################################################################
#                                                                         #
#             And write a nice configuration summary to file:             #
#                                                                         #
###########################################################################

SET(_log "${CMAKE_BINARY_DIR}/summary.log")

FILE(WRITE ${_log}
"###
#
#  ${DEAL_II_PACKAGE_NAME} configuration:
#
#        CMAKE_BUILD_TYPE:       ${CMAKE_BUILD_TYPE}
#        BUILD_SHARED_LIBS:      ${BUILD_SHARED_LIBS}
#        CMAKE_INSTALL_PREFIX:   ${CMAKE_INSTALL_PREFIX}
#        CMAKE_SOURCE_DIR:       ${CMAKE_SOURCE_DIR} (Version ${DEAL_II_PACKAGE_VERSION})
#        CMAKE_BINARY_DIR:       ${CMAKE_BINARY_DIR}
#        CMAKE_CXX_COMPILER:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME}
#                                ${CMAKE_CXX_COMPILER}
#
#  Compiler flags used for this build:
#        CMAKE_CXX_FLAGS:                     ${CMAKE_CXX_FLAGS}
"
  )
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  FILE(APPEND ${_log} "#        DEAL_II_CXX_FLAGS_RELEASE:           ${DEAL_II_CXX_FLAGS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  FILE(APPEND ${_log} "#        DEAL_II_CXX_FLAGS_DEBUG:             ${DEAL_II_CXX_FLAGS_DEBUG}\n")
ENDIF()
FILE(APPEND ${_log} "#        CMAKE_SHARED_LINKER_FLAGS:           ${CMAKE_SHARED_LINKER_FLAGS}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  FILE(APPEND ${_log} "#        DEAL_II_SHARED_LINKER_FLAGS_RELEASE: ${DEAL_II_SHARED_LINKER_FLAGS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  FILE(APPEND ${_log} "#        DEAL_II_SHARED_LINKER_FLAGS_DEBUG:   ${DEAL_II_SHARED_LINKER_FLAGS_DEBUG}\n")
ENDIF()
FILE(APPEND ${_log} "#        DEAL_II_SHARED_LINKER_FLAGS_DEBUG:   ${DEAL_II_SHARED_LINKER_FLAGS_DEBUG}\n")

IF(NOT DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
  FILE(APPEND ${_log} "#\n#  WARNING: DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS is set to OFF\n")
ELSEIF(NOT DEAL_II_KNOWN_COMPILER)
  FILE(APPEND ${_log} "#\n#  WARNING: Unknown compiler! Please set compiler flags by hand.\n")
ENDIF()

FILE(APPEND ${_log}
  "#\n#  Configured Features ("
  )
IF(DEFINED DEAL_II_ALLOW_BUNDLED)
  FILE(APPEND ${_log}
    "DEAL_II_ALLOW_BUNDLED = ${DEAL_II_ALLOW_BUNDLED}, "
    )
ENDIF()
IF(DEAL_II_FORCE_AUTODETECTION)
  FILE(APPEND ${_log}
    "!!! DEAL_II_FORCE_AUTODETECTION=ON !!!, "
    )
ENDIF()
FILE(APPEND ${_log}
  "DEAL_II_ALLOW_AUTODETECTION = ${DEAL_II_ALLOW_AUTODETECTION}):\n"
  )

GET_CMAKE_PROPERTY(_res VARIABLES)
FOREACH(_var ${_res})
  IF(_var MATCHES "DEAL_II_WITH")
    IF(${${_var}})
      # FEATURE is enabled
      STRING(REGEX REPLACE "^DEAL_II_WITH_" "" _feature ${_var})
      IF(FEATURE_${_feature}_EXTERNAL_CONFIGURED)
        FILE(APPEND ${_log} "#        ${_var} set up with external dependencies\n")
      ELSEIF(FEATURE_${_feature}_BUNDLED_CONFIGURED)
        IF(DEAL_II_FORCE_BUNDLED_${_feature})
          FILE(APPEND ${_log} "#        ${_var} set up with bundled packages (forced)\n")
        ELSE()
          FILE(APPEND ${_log} "#        ${_var} set up with bundled packages\n")
        ENDIF()
      ENDIF()
    ELSE()
      # FEATURE is disabled
      FILE(APPEND ${_log} "#      ( ${_var} = ${${_var}} )\n")
    ENDIF()
  ENDIF()
ENDFOREACH()

FILE(APPEND ${_log}
  "#\n#  Component configuration:\n"
  )
FOREACH(_var ${_res})
  IF(_var MATCHES "DEAL_II_COMPONENT")
    IF(${${_var}})
      FILE(APPEND ${_log} "#        ${_var}\n")
      STRING(REPLACE "DEAL_II_COMPONENT_" "" _component ${_var})
      LIST(APPEND _components ${_component})
    ELSE()
      FILE(APPEND ${_log} "#      ( ${_var} = ${${_var}} )\n")
    ENDIF()
  ENDIF()
ENDFOREACH()

FILE(APPEND ${_log}
  "#\n###"
  )
