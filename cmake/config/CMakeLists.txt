## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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

#
# This file sets up the project configuration consisting of
#
#   deal.IIConfig.cmake
#   deal.IIVersionConfig.cmake
#
# We support two configurations out of which deal.II can be used - directly
# from the build directory or after installation. So we have to prepare
# two distinct setups.
#

MESSAGE(STATUS "Setting up project configuration")

#
# Configure the template-arguments file
#
CONFIGURE_FILE( # for binary dir:
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/template-arguments
  )


########################################################################
#                                                                      #
#                   Setup and install cmake macros:                    #
#                                                                      #
########################################################################

SET(_macros
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_invoke_autopilot.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_target.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_query_git_information.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_add_test.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_pickup_tests.cmake
  )
FILE(COPY ${_macros}
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/macros
  )
INSTALL(FILES ${_macros}
  DESTINATION ${DEAL_II_SHARE_RELDIR}/macros
  COMPONENT library
  )


########################################################################
#                                                                      #
#        Configure and install the cmake project configuration:        #
#                                                                      #
########################################################################

#
# Do not force --as-needed for executables on user side:
#

STRIP_FLAG(DEAL_II_LINKER_FLAGS "-Wl,--as-needed")

#
# Strip -Wno-deprecated-declarations from DEAL_II_CXX_FLAGS so that
# deprecation warnings are actually shown for user code:
#
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wno-deprecated-declarations")

#
# Populate a bunch of CONFIG_* variables with useful information:
#

FOREACH(_build ${DEAL_II_BUILD_TYPES})

  IF(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    SET(_keyword "general")
  ELSE()
    IF(_build MATCHES DEBUG)
      SET(_keyword "debug")
    ELSE()
      SET(_keyword "optimized")
    ENDIF()
  ENDIF()

  #
  # Build up library name depending on link type and platform:
  #

  IF(BUILD_SHARED_LIBS)
    SET(_type "SHARED")
  ELSE()
    SET(_type "STATIC")
  ENDIF()

  SET(_name "${CMAKE_${_type}_LIBRARY_PREFIX}${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}${CMAKE_${_type}_LIBRARY_SUFFIX}")
  SET(CONFIG_LIBRARIES_${_build}
    "\${DEAL_II_PATH}/${DEAL_II_LIBRARY_RELDIR}/${_name}"
    ${DEAL_II_LIBRARIES_${_build}}
    ${DEAL_II_LIBRARIES}
    )
  LIST(APPEND CONFIG_LIBRARIES ${_keyword} \${DEAL_II_LIBRARIES_${_build}})

  SET(CONFIG_TARGET_${_build} ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX})
  LIST(APPEND CONFIG_TARGET ${_keyword} \${DEAL_II_TARGET_${_build}})
ENDFOREACH()

#
# For binary dir:
#

SET(CONFIG_BUILD_DIR TRUE)
SET(CONFIG_INCLUDE_DIRS
  \${DEAL_II_PATH}/include
  ${CMAKE_SOURCE_DIR}/include/
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_USER_INCLUDE_DIRS}
  )
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )

#
# For installation:
#

SET(CONFIG_BUILD_DIR FALSE)
SET(CONFIG_INCLUDE_DIRS
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
  ${DEAL_II_USER_INCLUDE_DIRS}
  )
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )
INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

#
# Append feature configuration to all configuration files:
#

SET(_files
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  )
FOREACH(_file ${_files})
  FILE(APPEND ${_file} "\n\n#\n# Feature configuration:\n#\n\n")
ENDFOREACH()

GET_CMAKE_PROPERTY(_res VARIABLES)
FOREACH(_var ${_res})
  IF(_var MATCHES "DEAL_II_.*WITH")
    LIST(APPEND _additional_config_variables ${_var})
  ENDIF()
ENDFOREACH()

SET(_deal_ii_features_sorted ${DEAL_II_FEATURES})
LIST(SORT _deal_ii_features_sorted)
FOREACH(_name ${_deal_ii_features_sorted})
  SET(_var DEAL_II_WITH_${_name})

  FOREACH(_file ${_files})
    FILE(APPEND ${_file} "SET(${_var} ${${_var}})\n")
    IF(${_var} AND NOT "${${_name}_VERSION}" STREQUAL "")
      FILE(APPEND ${_file}
        "SET(DEAL_II_${_name}_VERSION \"${${_name}_VERSION}\")\n"
        )
    ENDIF()
    FOREACH(_additional ${_additional_config_variables})
      IF(_additional MATCHES "DEAL_II_${_name}_WITH")
        FILE(APPEND ${_file} "SET(${_additional} ${${_additional}})\n")
      ENDIF()
    ENDFOREACH()
  ENDFOREACH()
ENDFOREACH()

MESSAGE(STATUS "Setting up project configuration - Done")
