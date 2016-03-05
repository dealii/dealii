## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2016 by the deal.II authors
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
#                Query for git repository information:                 #
#                                                                      #
########################################################################

DEAL_II_QUERY_GIT_INFORMATION()

FILE(WRITE ${CMAKE_BINARY_DIR}/revision.log
"###
#
#  Git information:
#        Branch:   ${DEAL_II_GIT_BRANCH}
#        Revision: ${DEAL_II_GIT_REVISION}
#
###"
  )


########################################################################
#                                                                      #
#              Write a nice configuration summary to file:             #
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
#        CMAKE_SOURCE_DIR:       ${CMAKE_SOURCE_DIR}
"
  )
IF("${DEAL_II_GIT_SHORTREV}" STREQUAL "")
  _both("#                                (version ${DEAL_II_PACKAGE_VERSION})\n")
ELSE()
  _both("#                                (version ${DEAL_II_PACKAGE_VERSION}, shortrev ${DEAL_II_GIT_SHORTREV})\n")
ENDIF()
_both(
"#        CMAKE_BINARY_DIR:       ${CMAKE_BINARY_DIR}
#        CMAKE_CXX_COMPILER:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}
#                                ${CMAKE_CXX_COMPILER}
"
  )

IF(CMAKE_C_COMPILER_WORKS)
  _detailed("#        CMAKE_C_COMPILER:       ${CMAKE_C_COMPILER}\n")
ENDIF()
IF(CMAKE_Fortran_COMPILER_WORKS)
  _detailed("#        CMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}\n")
ENDIF()
_detailed("#        CMAKE_GENERATOR:        ${CMAKE_GENERATOR}\n")

IF(CMAKE_CROSSCOMPILING)
  _both(
    "#\n#        CROSSCOMPILING!\n"
    )
ENDIF()

IF(DEAL_II_STATIC_EXECUTABLE)
  _both(
    "#\n#        STATIC LINKAGE!\n"
    )
ENDIF()

_both("#\n")

_detailed(
"#  Base configuration (prior to feature configuration):
#        DEAL_II_CXX_FLAGS:            ${BASE_CXX_FLAGS}
"
  )
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_CXX_FLAGS_RELEASE:    ${BASE_CXX_FLAGS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_CXX_FLAGS_DEBUG:      ${BASE_CXX_FLAGS_DEBUG}\n")
ENDIF()

_detailed("#        DEAL_II_LINKER_FLAGS:         ${BASE_LINKER_FLAGS}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_LINKER_FLAGS_RELEASE: ${BASE_LINKER_FLAGS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_LINKER_FLAGS_DEBUG:   ${BASE_LINKER_FLAGS_DEBUG}\n")
ENDIF()

_detailed("#        DEAL_II_DEFINITIONS:          ${BASE_DEFINITIONS}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_DEFINITIONS_RELEASE:  ${BASE_DEFINITIONS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_DEFINITIONS_DEBUG:    ${BASE_DEFINITIONS_DEBUG}\n")
ENDIF()

_detailed("#        DEAL_II_USER_DEFINITIONS:     ${BASE_DEFINITIONS}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_USER_DEFINITIONS_REL: ${BASE_DEFINITIONS_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_USER_DEFINITIONS_DEB: ${BASE_DEFINITIONS_DEBUG}\n")
ENDIF()

_detailed("#        DEAL_II_INCLUDE_DIRS          ${BASE_INCLUDE_DIRS}\n")
_detailed("#        DEAL_II_USER_INCLUDE_DIRS:    ${BASE_USER_INCLUDE_DIRS}\n")
_detailed("#        DEAL_II_BUNDLED_INCLUDE_DIRS: ${BASE_BUNDLED_INCLUDE_DIRS}\n")

_detailed("#        DEAL_II_LIBRARIES:            ${BASE_LIBRARIES}\n")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_LIBRARIES_RELEASE:    ${BASE_LIBRARIES_RELEASE}\n")
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_LIBRARIES_DEBUG:      ${BASE_LIBRARIES_DEBUG}\n")
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


SET(_deal_ii_features_sorted ${DEAL_II_FEATURES})
LIST(SORT _deal_ii_features_sorted)
FOREACH(_feature ${_deal_ii_features_sorted})
  SET(_var DEAL_II_WITH_${_feature})

  IF(${${_var}})
    #
    # The feature is enabled:
    #
    IF(FEATURE_${_feature}_EXTERNAL_CONFIGURED)
      _both("#        ${_var} set up with external dependencies\n")
    ELSEIF(FEATURE_${_feature}_BUNDLED_CONFIGURED)
      IF(DEAL_II_FORCE_BUNDLED_${_feature})
        _both("#        ${_var} set up with bundled packages (forced)\n")
      ELSE()
        _both("#        ${_var} set up with bundled packages\n")
      ENDIF()
    ELSE()
     _both("#        ${_var} = ${${_var}}\n")
    ENDIF()

    #
    # Print out version number:
    #
    IF(DEFINED ${_feature}_VERSION)
      _detailed("#            ${_feature}_VERSION = ${${_feature}_VERSION}\n")
    ENDIF()

    #
    # Special version numbers:
    #
    IF(_feature MATCHES "THREADS" AND DEFINED TBB_VERSION)
      _detailed("#            TBB_VERSION = ${TBB_VERSION}\n")
    ENDIF()
    IF(_feature MATCHES "MPI" AND DEFINED OMPI_VERSION)
      _detailed("#            OMPI_VERSION = ${OMPI_VERSION}\n")
    ENDIF()

    #
    # Print out ${_feature}_DIR:
    #
    IF(NOT "${${_feature}_DIR}" STREQUAL "")
      _detailed("#            ${_feature}_DIR = ${${_feature}_DIR}\n")
    ENDIF()

    #
    # Print the feature configuration:
    #
    FOREACH(_var2
      C_COMPILER CXX_COMPILER Fortran_COMPILER
      ${DEAL_II_STRING_SUFFIXES} ${DEAL_II_LIST_SUFFIXES}
      )
      IF(DEFINED ${_feature}_${_var2})
        _detailed("#            ${_feature}_${_var2} = ${${_feature}_${_var2}}\n")
      ENDIF()
    ENDFOREACH()
  ELSE()
    # FEATURE is disabled
    _both("#      ( ${_var} = ${${_var}} )\n")
  ENDIF()
ENDFOREACH()

_both(
  "#\n#  Component configuration:\n"
  )

FOREACH(_component ${DEAL_II_COMPONENTS})
  SET(_var DEAL_II_COMPONENT_${_component})
  IF(${${_var}})
    _both("#        ${_var}\n")
  ELSE()
    _both("#      ( ${_var} = ${${_var}} )\n")
  ENDIF()
ENDFOREACH()

_summary(
"#\n#  Detailed information (compiler flags, feature configuration) can be found in detailed.log
#\n#  Run  $ "
  )
IF(CMAKE_GENERATOR MATCHES "Ninja")
  _summary("ninja ")
ELSE()
_summary("make ")
ENDIF()
_summary("info  to print a help message with a list of top level targets\n")

_both("#\n###")
