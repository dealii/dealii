## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# This file sets up the project installation and configuration consisting of
#
#   deal.IIConfig.cmake
#   deal.IIVersionConfig.cmake
#
# We support two configurations out of which deal.II can be used - directly
# from the build directory or after installation. So we have to prepare
# two distinct setups.
#


message(STATUS "Setting up library installation")

set(_targets ${DEAL_II_TARGET_NAME})
foreach(build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER ${build} build_lowercase)
  list(APPEND _targets ${DEAL_II_TARGET_NAME}_${build_lowercase})
  if(DEAL_II_WITH_CXX20_MODULE)
    list(APPEND _targets ${DEAL_II_TARGET_NAME}_module_${build_lowercase})
  endif()
endforeach()

if(CMAKE_VERSION VERSION_LESS 3.28)
  install(TARGETS ${_targets} ${DEAL_II_TARGETS} ${DEAL_II_TARGETS_DEBUG} ${DEAL_II_TARGETS_RELEASE}
    COMPONENT library
    EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
    RUNTIME DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
    LIBRARY DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    ARCHIVE DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    )
else()
  install(TARGETS ${_targets} ${DEAL_II_TARGETS} ${DEAL_II_TARGETS_DEBUG} ${DEAL_II_TARGETS_RELEASE}
    COMPONENT library
    EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
    RUNTIME DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
    LIBRARY DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    ARCHIVE DESTINATION ${DEAL_II_LIBRARY_RELDIR}
    FILE_SET external_module_interface_units
      DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/module
    FILE_SET module_interface_units
      DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/module
    FILE_SET headers
      DESTINATION ${DEAL_II_INCLUDE_RELDIR}
    )
endif()

install(EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
  NAMESPACE "${DEAL_II_TARGET_NAME}::"
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR})
export(TARGETS ${_targets} ${DEAL_II_TARGETS} ${DEAL_II_TARGETS_DEBUG} ${DEAL_II_TARGETS_RELEASE}
  NAMESPACE "${DEAL_II_TARGET_NAME}::"
  FILE ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
)

message(STATUS "Setting up library installation - Done")


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
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_add_test.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_invoke_autopilot.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_pickup_tests.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_query_git_information.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_target.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_shell_escape_option_groups.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_target_compile_flags.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_target_link_flags.cmake
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
# Populate a bunch of CONFIG_* variables with useful information:
#

foreach(_build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER ${_build} _build_lowercase)

  set(_keyword)
  if("${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    if(_build MATCHES DEBUG)
      set(_keyword "debug")
    else()
      set(_keyword "optimized")
    endif()
  endif()

  set(CONFIG_TARGET_${_build} ${DEAL_II_TARGET_NAME}::${DEAL_II_TARGET_NAME}_${_build_lowercase})
  list(APPEND CONFIG_TARGET ${_keyword} \${DEAL_II_TARGET_${_build}})
endforeach()

#
# Collect feature configuration that will be added to the configuration files:
#

set(FEATURE_CONFIGURATION_STRING "
#
# Feature configuration:
#

set(DEAL_II_WITH_CXX11 ON)
set(DEAL_II_WITH_CXX14 ON)
set(DEAL_II_WITH_CXX17 ON)
set(DEAL_II_WITH_CXX20 ${DEAL_II_HAVE_CXX20})
set(DEAL_II_WITH_THREADS ON)\n"
  )

set(_additional_config_variables)
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

  string(APPEND FEATURE_CONFIGURATION_STRING "set(${_var} ${${_var}})\n")
  #
  # Do not pollute deal.IIConfig.cmake with package details of
  # unconfigured features.
  #
  if(${_var})
    if(NOT "${${_name}_VERSION}" STREQUAL "")
      string(APPEND FEATURE_CONFIGURATION_STRING
        "set(DEAL_II_${_name}_VERSION \"${${_name}_VERSION}\")\n"
        )
    endif()
    foreach(_additional ${_additional_config_variables})
      if(_additional MATCHES "(DEAL_II_${_name}_WITH|DEAL_II_FEATURE_${_name})")
        string(APPEND FEATURE_CONFIGURATION_STRING
          "set(${_additional} ${${_additional}})\n"
          )
      endif()
    endforeach()
  endif()
endforeach()

#
# For binary dir:
#

set(CONFIG_BUILD_DIR TRUE)

#
# Filter out -Wno-deprecated-decarations so that user projects actually see
# deprecation warnings:
#
string(REPLACE "-Wno-deprecated-declarations" "" DEAL_II_WARNING_FLAGS "${DEAL_II_WARNING_FLAGS}")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigGit.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigGit.cmake
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

set(CONFIG_BUILD_DIR FALSE)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigGit.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigGit.cmake
  @ONLY
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigGit.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

#
# Job's done.
#

message(STATUS "Setting up project configuration - Done")
