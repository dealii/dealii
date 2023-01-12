## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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

message(STATUS "Setting up project configuration")

#
# Configure the template-arguments file
#
configure_file( # for binary dir:
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/template-arguments
  )


########################################################################
#                                                                      #
#                   Setup and install cmake macros:                    #
#                                                                      #
########################################################################

set(_macros
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_invoke_autopilot.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_target.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_query_git_information.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_add_test.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_pickup_tests.cmake
  )
file(COPY ${_macros}
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/macros
  )
install(FILES ${_macros}
  DESTINATION ${DEAL_II_SHARE_RELDIR}/macros
  COMPONENT library
  )


########################################################################
#                                                                      #
#        Configure and install the cmake project configuration:        #
#                                                                      #
########################################################################

#
# Strip -Wno-deprecated-declarations from DEAL_II_CXX_FLAGS so that
# deprecation warnings are actually shown for user code:
#

strip_flag(DEAL_II_CXX_FLAGS "-Wno-deprecated-declarations")

#
# Strip -Wsuggest-override. Enabling this warning is a stylistic choice
# that might create a substantial amount of warnings in user codebases.
# Let's be nice and disable the flag for now.
#

strip_flag(DEAL_II_CXX_FLAGS "-Wsuggest-override")

#
# Populate a bunch of CONFIG_* variables with useful information:
#

foreach(_build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER ${_build} _build_lowercase)

  if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    set(_keyword "general")
  else()
    if(_build MATCHES DEBUG)
      set(_keyword "debug")
    else()
      set(_keyword "optimized")
    endif()
  endif()

  #
  # Build up library name depending on link type and platform:
  #

  if(BUILD_SHARED_LIBS)
    set(_type "SHARED")
  else()
    set(_type "STATIC")
  endif()

  set(_name "${CMAKE_${_type}_LIBRARY_PREFIX}${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}${CMAKE_${_type}_LIBRARY_SUFFIX}")
  set(CONFIG_LIBRARIES_${_build}
    "\${DEAL_II_PATH}/${DEAL_II_LIBRARY_RELDIR}/${_name}"
    ${DEAL_II_LIBRARIES_${_build}}
    ${DEAL_II_LIBRARIES}
    )
  list(APPEND CONFIG_LIBRARIES ${_keyword} \${DEAL_II_LIBRARIES_${_build}})

  set(CONFIG_TARGET_${_build} ${DEAL_II_NAMESPACE}::${DEAL_II_NAMESPACE}_${_build_lowercase})
  list(APPEND CONFIG_TARGET ${_keyword} \${DEAL_II_TARGET_${_build}})
endforeach()

#
# For binary dir:
#

set(CONFIG_BUILD_DIR TRUE)
set(CONFIG_INCLUDE_DIRS
  \${DEAL_II_PATH}/include
  ${CMAKE_SOURCE_DIR}/include/
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_INCLUDE_DIRS}
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )

#
# For installation:
#

# Make sure that we only ever record the bundled subdirectory if it is in
# fact created:
if("${DEAL_II_BUNDLED_INCLUDE_DIRS}" STREQUAL "")
  set(_installed_bundled "")
else()
  set(_installed_bundled
    "\${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled"
    )
endif()

set(CONFIG_BUILD_DIR FALSE)
set(CONFIG_INCLUDE_DIRS
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}
  ${_installed_bundled}
  ${DEAL_II_INCLUDE_DIRS}
  )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

#
# Append feature configuration to all configuration files:
#

set(_files
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  )
foreach(_file ${_files})
  file(APPEND ${_file} "\n\n#\n# Feature configuration:\n#\n\n")
  file(APPEND ${_file}
    "set(DEAL_II_WITH_CXX11 ON)\n"
    "set(DEAL_II_WITH_CXX14 ON)\n"
    "set(DEAL_II_WITH_CXX17 ${DEAL_II_HAVE_CXX17})\n"
    "set(DEAL_II_WITH_THREADS ON)\n"
    )
endforeach()

get_cmake_property(_res VARIABLES)
foreach(_var ${_res})
  if(_var MATCHES "DEAL_II_.*(WITH|FEATURE.*BUNDLED)")
    list(APPEND _additional_config_variables ${_var})
  endif()
endforeach()

set(_deal_ii_features_sorted ${DEAL_II_FEATURES})
list(SORT _deal_ii_features_sorted)
foreach(_name ${_deal_ii_features_sorted})
  set(_var DEAL_II_WITH_${_name})

  foreach(_file ${_files})
    file(APPEND ${_file} "set(${_var} ${${_var}})\n")
    #
    # Do not pollute deal.IIConfig.cmake with package details of
    # unconfigured features.
    #
    if(${_var})
      if(NOT "${${_name}_VERSION}" STREQUAL "")
        file(APPEND ${_file}
          "set(DEAL_II_${_name}_VERSION \"${${_name}_VERSION}\")\n"
          )
      endif()

      foreach(_additional ${_additional_config_variables})
        if(_additional MATCHES "(DEAL_II_${_name}_WITH|DEAL_II_FEATURE_${_name})")
          file(APPEND ${_file} "set(${_additional} ${${_additional}})\n")
        endif()
      endforeach()
    endif()
  endforeach()
endforeach()

#
# Job's done.
#

message(STATUS "Setting up project configuration - Done")

