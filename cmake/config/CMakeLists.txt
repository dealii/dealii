## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2019 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
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
# Strip -Wsuggest-override. Enabling this warning is a stylistic choice
# that might create a substantial amount of warnings in user codebases.
# Let's be nice and disable the flag for now.
#

STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wsuggest-override")

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

# Make sure that we only ever record the bundled subdirectory if it is in
# fact created:
IF("${DEAL_II_BUNDLED_INCLUDE_DIRS}" STREQUAL "")
  SET(_installed_bundled "")
ELSE()
  SET(_installed_bundled
    "\${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled"
    )
ENDIF()

SET(CONFIG_BUILD_DIR FALSE)
SET(CONFIG_INCLUDE_DIRS
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}
  ${_installed_bundled}
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
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}FeatureConfig.cmake
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
  FILE(APPEND ${_file}
    "SET(DEAL_II_WITH_CXX11 ON)\n"
    "SET(DEAL_II_WITH_CXX14 ON)\n"
    "SET(DEAL_II_WITH_CXX17 ${DEAL_II_HAVE_CXX17})\n"
    "SET(DEAL_II_WITH_THREADS ON)\n"
    )
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
    #
    # Do not pollute deal.IIConfig.cmake with package details of
    # unconfigured features.
    #
    IF(${_var})
      IF(NOT "${${_name}_VERSION}" STREQUAL "")
        FILE(APPEND ${_file}
          "SET(DEAL_II_${_name}_VERSION \"${${_name}_VERSION}\")\n"
          )
      ENDIF()
      FOREACH(_additional ${_additional_config_variables})
        IF(_additional MATCHES "DEAL_II_${_name}_WITH")
          FILE(APPEND ${_file} "SET(${_additional} ${${_additional}})\n")
        ENDIF()
      ENDFOREACH()
    ENDIF()
  ENDFOREACH()
ENDFOREACH()


########################################################################
#                                                                      #
#                           pkgconfig files                            #
#                                                                      #
########################################################################

#
# Set up the pkgconfig configuration files consisting of
#
#   deal.II_debug      # for the debug variant of the library
#   deal.II_release    # for the release variant of the library
#
# Similarly to the CMake project configuration, we provide pkgconfig files
# directly for the build directory, as well as for the final installation.
# So we have to prepare two distinct setups.
#
# pkgconfig looks for *.pc files in an environmental variable called
# PKG_CONFIG_PATH. So, to use the library in the build directory issue
# export PKG_CONFIG_PATH=/path/to/BUILD_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
#
# To use the library in the installed location
# export PKG_CONFIG_PATH=/path/to/INSTALL_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
#

#
# We need to gather some variables for the regex below to work.
#

SET(_library_prefixes "")
IF (CMAKE_SHARED_LIBRARY_PREFIX)
  LIST(APPEND _library_prefixes ${CMAKE_SHARED_LIBRARY_PREFIX})
ENDIF()
IF (CMAKE_STATIC_LIBRARY_PREFIX)
  LIST(APPEND _library_prefixes ${CMAKE_STATIC_LIBRARY_PREFIX})
ENDIF()
STRING(REPLACE ";" "|" _library_prefixes "${_library_prefixes}")

SET(_library_suffixes "")
IF (CMAKE_SHARED_LIBRARY_SUFFIX)
  LIST(APPEND _library_suffixes ${CMAKE_SHARED_LIBRARY_SUFFIX})
ENDIF()
IF (CMAKE_STATIC_LIBRARY_SUFFIX)
  LIST(APPEND _library_suffixes ${CMAKE_STATIC_LIBRARY_SUFFIX})
ENDIF()
STRING(REPLACE ";" "|" _library_suffixes "${_library_suffixes}")

#
# Build up the link line from our list of libraries:
#

FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER ${_build} _build_lowercase)

  SET(_name "${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}")

  SET(CONFIG_RPATH_${_build}_PC "\\\${libdir}")
  SET(CONFIG_LIBRARIES_${_build}_PC "-L\${libdir} -l${_name}")

  FOREACH(_lib ${DEAL_II_LIBRARIES_${_build}} ${DEAL_II_LIBRARIES})

    GET_FILENAME_COMPONENT(_name ${_lib} NAME)
    GET_FILENAME_COMPONENT(_dir ${_lib} PATH)

    IF("${_dir}" STREQUAL "")
      # ${_lib} is a simple library name, just add it to the link line:
      SET(_library_string "-l${_lib}")

    ELSE()
      # ${_lib} is a full library path:

      #
      # Only append a library directory if it is new ...
      #

      LIST(FIND CONFIG_RPATH_${_build}_PC "${_dir}" _index)
      IF (${_index} EQUAL -1)
        SET(_library_string "-L${_dir} ")
        LIST(APPEND CONFIG_RPATH_${_build}_PC ${_dir})
      ELSE()
        SET(_library_string "")
      ENDIF()

      # Recover short name:
      IF(_library_prefixes)
        STRING(REGEX REPLACE "^(${_library_prefixes})" "" _name "${_name}")
      ENDIF()
      IF(_library_suffixes)
        STRING(REGEX REPLACE "(${_library_suffixes})$" "" _name "${_name}")
      ENDIF()
      SET(_library_string "${_library_string}-l${_name}")
    ENDIF()

    SET(CONFIG_LIBRARIES_${_build}_PC
      "${CONFIG_LIBRARIES_${_build}_PC} ${_library_string}"
      )
  ENDFOREACH()

  TO_STRING_AND_ADD_PREFIX(CONFIG_RPATH_${_build}_PC
    "-Wl,-rpath," ${CONFIG_RPATH_${_build}_PC}
    )
ENDFOREACH()


#
# For binary dir (aka build dir):
#

SET(CONFIG_PATH_PC "${CMAKE_BINARY_DIR}")

TO_STRING_AND_ADD_PREFIX(CONFIG_INCLUDE_DIRS_PC "-I"
  \\\${prefix}/include
  ${CMAKE_SOURCE_DIR}/include
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_USER_INCLUDE_DIRS}
  )

FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER ${_build} _build_lowercase)

  SET(_config_directory "${CMAKE_BINARY_DIR}/${DEAL_II_PKGCONF_RELDIR}")

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/config_${_build_lowercase}.pc.in
    ${_config_directory}/${DEAL_II_PROJECT_CONFIG_NAME}_${_build_lowercase}.pc
    @ONLY
    )
ENDFOREACH()

#
# For installation:
#

SET(CONFIG_PATH_PC "${CMAKE_INSTALL_PREFIX}")

TO_STRING_AND_ADD_PREFIX(CONFIG_INCLUDE_DIRS_PC "-I"
  \\\${includedir}
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_USER_INCLUDE_DIRS}
  )

FOREACH(_build ${DEAL_II_BUILD_TYPES})
  STRING(TOLOWER ${_build} _build_lowercase)

  #
  # Only populate the pkgconf files for the install directory if
  # CMAKE_INSTALL_RPATH_USE_LINK_PATH is true.
  #
  # We use this a heuristic for now to decide whether the user actually
  # wants to have RPATHs in configuration files after installation.
  #
  # FIXME: Unify RPATH handling between cmake and pkgconf configuration and
  # clearly document RPATH behavior.
  #
  IF(NOT CMAKE_INSTALL_RPATH_USE_LINK_PATH)
    SET(CONFIG_RPATH_${_build}_PC "")
  ENDIF()

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/config_${_build_lowercase}.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_${_build_lowercase}.pc
    @ONLY
    )
  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_${_build_lowercase}.pc
    DESTINATION ${DEAL_II_PKGCONF_RELDIR}
    COMPONENT library
    )
ENDFOREACH()

#
# Job's done.
#

MESSAGE(STATUS "Setting up project configuration - Done")

