## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------


########################################################################
#                                                                      #
#                      Finalize the configuration:                     #
#                                                                      #
########################################################################

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
# Deduplicate entries in DEAL_II_USER_INCLUDE_DIRS
#
IF(NOT "${DEAL_II_USER_INCLUDE_DIRS}" STREQUAL "")
  LIST(REMOVE_DUPLICATES DEAL_II_USER_INCLUDE_DIRS)
ENDIF()

#
# Deduplicate entries in DEAL_II_EXTERNAL_LIBRARIES(_...)
# in reverse order:
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

#
# Cleanup some files used for storing the names of all object targets that
# will be bundled to the deal.II library.
# (Right now, i.e. cmake 2.8.8, this is the only reliable way to get
# information into a global scope...)
#
FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER "${_build}" _build_lowercase)
  FILE(REMOVE
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_objects_${_build_lowercase}
    )
ENDFOREACH()
FILE(WRITE
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_source_includes
  ""
  )

#
# Cleanup deal.IITargets.cmake in the build directory:
#
FILE(REMOVE
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}BuildTargets.cmake
  )


########################################################################
#                                                                      #
#            And write a nice configuration summary to file:           #
#                                                                      #
########################################################################

SET(_log_detailed "${CMAKE_BINARY_DIR}/detailed.log")
SET(_log_summary  "${CMAKE_BINARY_DIR}/summary.log")
FILE(REMOVE ${_log_detailed} ${_log_summary})

MACRO(_both)
  # Write to both log files:
  FILE(APPEND ${_log_detailed} "${ARGN}")
  FILE(APPEND ${_log_summary} "${ARGN}")
ENDMACRO()
MACRO(_detailed)
  # Only write to detailed.log:
  FILE(APPEND ${_log_detailed} "${ARGN}")
ENDMACRO()
MACRO(_summary)
  # Only write to summary.log:
  FILE(APPEND ${_log_summary} "${ARGN}")
ENDMACRO()


_both(
"###
#
#  ${DEAL_II_PACKAGE_NAME} configuration:
#        CMAKE_BUILD_TYPE:       ${CMAKE_BUILD_TYPE}
#        BUILD_SHARED_LIBS:      ${BUILD_SHARED_LIBS}
#        CMAKE_INSTALL_PREFIX:   ${CMAKE_INSTALL_PREFIX}
#        CMAKE_SOURCE_DIR:       ${CMAKE_SOURCE_DIR} (Version ${DEAL_II_PACKAGE_VERSION})
#        CMAKE_BINARY_DIR:       ${CMAKE_BINARY_DIR}
#        CMAKE_CXX_COMPILER:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}
#                                ${CMAKE_CXX_COMPILER}
"
  )
IF(CMAKE_CROSSCOMPILING)
  _both(
    "#\n#        CROSSCOMPILING!\n"
    "#        DEAL_II_NATIVE:         ${DEAL_II_NATIVE}\n"
    )
ENDIF()

IF(DEAL_II_STATIC_EXECUTABLE)
  _both(
    "#\n#        STATIC LINKAGE!\n"
    )
ENDIF()

_both("#\n")

_detailed(
"#  Compiler flags used for this build:
#        CMAKE_CXX_FLAGS:              ${CMAKE_CXX_FLAGS}
"
  )
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_CXX_FLAGS_RELEASE:    ${DEAL_II_CXX_FLAGS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_CXX_FLAGS_DEBUG:      ${DEAL_II_CXX_FLAGS_DEBUG}\n")
ENDIF()
_detailed("#        DEAL_II_LINKER_FLAGS:         ${DEAL_II_LINKER_FLAGS}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_LINKER_FLAGS_RELEASE: ${DEAL_II_LINKER_FLAGS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_LINKER_FLAGS_DEBUG:   ${DEAL_II_LINKER_FLAGS_DEBUG}\n")
ENDIF()
_detailed("#\n")

IF(NOT DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS)
  _both("#  WARNING: DEAL_II_SETUP_DEFAULT_COMPILER_FLAGS is set to OFF\n")
ENDIF()
_both("#  Configured Features (")
IF(DEFINED DEAL_II_ALLOW_BUNDLED)
  _both("DEAL_II_ALLOW_BUNDLED = ${DEAL_II_ALLOW_BUNDLED}, ")
ENDIF()
IF(DEAL_II_FORCE_AUTODETECTION)
  _both("!!! DEAL_II_FORCE_AUTODETECTION=ON !!!, ")
ENDIF()
_both("DEAL_II_ALLOW_AUTODETECTION = ${DEAL_II_ALLOW_AUTODETECTION}):\n")


#
# Cache for quicker access to avoid the O(n^2) complexity of a loop over
# _all_ defined variables.
#
GET_CMAKE_PROPERTY(_variables VARIABLES)
FOREACH(_var ${_variables})
  IF(_var MATCHES "DEAL_II_WITH")
    LIST(APPEND _features "${_var}")
  ELSEIF(_var MATCHES "DEAL_II_COMPONENT")
    LIST(APPEND _components "${_var}")
  ELSEIF(_var MATCHES "(MPI_CXX_COMPILER|MPI_CXX_COMPILE_FLAGS|MPI_CXX_LINK_FLAGS)")
    LIST(APPEND _features_config ${_var})
  ELSEIF(_var MATCHES "(LIBRARIES|INCLUDE_PATH|INCLUDE_DIRS|LINKER_FLAGS)")
    LIST(APPEND _features_config ${_var})
  ENDIF()
ENDFOREACH()

FOREACH(_var ${_features})
  IF(${${_var}})
    # FEATURE is enabled
    STRING(REGEX REPLACE "^DEAL_II_WITH_" "" _feature ${_var})
    IF(FEATURE_${_feature}_EXTERNAL_CONFIGURED)
      _both("#        ${_var} set up with external dependencies\n")

      #
      # Print out version number:
      #
      IF(DEFINED ${_feature}_VERSION)
        _detailed("#            ${_feature}_VERSION = ${${_feature}_VERSION}\n")
      ENDIF()

      #
      # Print out ${_feature}_DIR:
      #
      IF(DEFINED ${_feature}_DIR)
        _detailed("#            ${_feature}_DIR = ${${_feature}_DIR}\n")
      ENDIF()

      #
      # Print the feature configuration:
      #
      FOREACH(_var2 ${_features_config})
        IF( # MPI:
            _var2 MATCHES "^${_feature}_CXX_(COMPILER|COMPILE_FLAGS|LINK_FLAGS|LIBRARIES|INCLUDE_PATH)$" OR
            # Boost:
            ( _feature MATCHES "BOOST" AND _var2 MATCHES "^Boost(_LIBRARIES|_INCLUDE_DIRS)$" ) OR
            # TBB:
            ( _feature MATCHES "THREADS" AND _var2 MATCHES "^TBB(_LIBRARIES|_INCLUDE_DIRS)$" ) OR
            # Generic:
            ( (NOT _var2 MATCHES "^(MPI|Boost)") AND
              _var2 MATCHES "^${_feature}_(INCLUDE_DIRS|LIBRARIES|LINKER_FLAGS)$" )
          )
          _detailed("#            ${_var2} = ${${_var2}}\n")
        ENDIF()
      ENDFOREACH()

    ELSEIF(FEATURE_${_feature}_BUNDLED_CONFIGURED)
      IF(DEAL_II_FORCE_BUNDLED_${_feature})
        _both("#        ${_var} set up with bundled packages (forced)\n")
      ELSE()
        _both("#        ${_var} set up with bundled packages\n")
      ENDIF()
    ELSE()
     _both("#        ${_var} = ${${_var}}\n")
    ENDIF()
  ELSE()
    # FEATURE is disabled
    _both("#      ( ${_var} = ${${_var}} )\n")
  ENDIF()
ENDFOREACH()

_both(
  "#\n#  Component configuration:\n"
  )
FOREACH(_var ${_components})
  IF(_var MATCHES "DEAL_II_COMPONENT")
    IF(${${_var}})
      _both("#        ${_var}\n")
      STRING(REPLACE "DEAL_II_COMPONENT_" "" _component ${_var})
      LIST(APPEND _components ${_component})
    ELSE()
      _both("#      ( ${_var} = ${${_var}} )\n")
    ENDIF()
  ENDIF()
ENDFOREACH()

_summary(
  "#\n# Detailed information (compiler flags, feature configuration) can be found in detailed.log\n"
  )

_both("#\n###")

