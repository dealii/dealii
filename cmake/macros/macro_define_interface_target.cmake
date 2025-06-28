## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2023 - 2025 by the deal.II authors
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
# define_interface_target(<feature> [target_name])
#
# This function defines interface targets for a given feature described by
# the CMake variables <FEATURE>_<SUFFIX>, where the suffix is one of the
# following:
#
#   LIBRARIES    ( LIBRARIES_RELEASE    LIBRARIES_DEBUG )
#   TARGETS      ( TARGETS_RELEASE      TARGETS_DEBUG )
#     - populating the LINK_LIBRARIES target property
#   INCLUDE_DIRS
#     - populating the INCLUDE_DIRECTORIES target property
#   DEFINITIONS  ( DEFINITIONS_RELEASE  DEFINITIONS_DEBUG )
#     - populating the COMPILE_DEFINITIONS target property
#   CXX_FLAGS    ( CXX_FLAGS_RELEASE    CXX_FLAGS_DEBUG )
#     - populating the COMPILE_OPTIONS target property
#   LINKER_FLAGS ( LINKER_FLAGS_RELEASE LINKER_FLAGS_DEBUG )
#     - populating the LINK_OPTIONS target property
#
# If any *_DEBUG or *_RELEASE variable is populated then this macro defines
# two targets named interface_<feature>_debug and
# interface_<feature>_release. If no *_DEBUG or *_RELEASE variants are
# defined then a single interface_<feature> target is defined.
#
#  - All interface targets are added automatically to
#      ${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
#    for build and install locations.
#
# - The function also adds all targets to DEAL_II_TARGETS (for non
#   debug/release split) and DEAL_II_TARGETS_(DEBUG|RELEASE) (for split
#   targets), respectively.
#
# The default name, i.e., interface_<feature>(|_debug|_release) can be
# overridden by the optional second argument. For example,
#   define_interface_target(DEAL_II base_configuration)
# will define interface_base_configuration* targets but query all
# information from DEAL_II_* variables.
#

function(define_interface_target _feature)
  string(TOLOWER "${_feature}" _feature_lowercase)
  if(NOT "${ARGN}" STREQUAL "")
    set(_feature_lowercase "${ARGN}")
  endif()

  set(_builds "dummy")
  if(${_feature}_SPLIT_CONFIGURATION)
    set(_builds ${DEAL_II_BUILD_TYPES})
  endif()

  foreach(_build ${_builds})
    set(_interface_target "interface_${_feature_lowercase}")

    if(${_feature}_SPLIT_CONFIGURATION)
      string(TOLOWER "${_build}" _build_lowercase)
      string(APPEND _interface_target "_${_build_lowercase}")
    endif()

    message(STATUS "")
    message(STATUS "Defining target: ${_interface_target}")

    add_library(${_interface_target} INTERFACE)

    if(DEFINED ${_feature}_VERSION)
      message(STATUS "    VERSION:             ${${_feature}_VERSION}")
      # CMake versions prior to 3.19 have a significantly more restrictive
      # set of allowed interface target properties.
      if(NOT CMAKE_VERSION VERSION_LESS 3.19)
        set_target_properties(${_interface_target}
          PROPERTIES VERSION "${${_feature_}_VERSION}"
          )
      endif()
    endif()

    set(_libraries)
    list(APPEND _libraries
      ${${_feature}_LIBRARIES} ${${_feature}_LIBRARIES_${_build}}
      )
    if(NOT "${_libraries}" STREQUAL "")
      foreach(_lib ${_libraries})
        #
        # Complain loudly if we encounter an undefined target:
        #
        if("${_lib}" MATCHES "::")
          message(FATAL_ERROR
            "Undefined imported target name \"${_lib}\" present when defining "
            "interface target \"${_interface_target}\""
            )
        endif()
      endforeach()

      message(STATUS "    LINK_LIBRARIES:      ${_libraries}")
      target_link_libraries(${_interface_target} INTERFACE ${_libraries})
    endif()

    if(DEFINED ${_feature}_INCLUDE_DIRS)
      message(STATUS "    INCLUDE_DIRECTORIES: ${${_feature}_INCLUDE_DIRS}")
      target_include_directories(${_interface_target}
        SYSTEM INTERFACE ${${_feature}_INCLUDE_DIRS}
        )
    endif()

    set(_definitions)
    list(APPEND _definitions ${${_feature}_DEFINITIONS} ${${_feature}_DEFINITIONS_${_build}})
    if(NOT "${_definitions}" STREQUAL "")
      message(STATUS "    COMPILE_DEFINITIONS: ${_definitions}")
      target_compile_definitions(${_interface_target} INTERFACE ${_definitions})
    endif()

    separate_arguments(_compile_options UNIX_COMMAND
      "${${_feature}_CXX_FLAGS} ${${_feature}_CXX_FLAGS_${_build}}"
      )
    shell_escape_option_groups(_compile_options)
    if(NOT "${_compile_options}" STREQUAL "")
      message(STATUS "    COMPILE_OPTIONS:     ${_compile_options}")
      target_compile_options(${_interface_target} INTERFACE ${_compile_options})
    endif()

    separate_arguments(_link_options UNIX_COMMAND
      "${${_feature}_LINKER_FLAGS} ${${_feature}_LINKER_FLAGS_${_build}}"
      )
    shell_escape_option_groups(_link_options)
    if(NOT "${_link_options}" STREQUAL "")
      message(STATUS "    LINK_OPTIONS:        ${_link_options}")
      target_link_options(${_interface_target} INTERFACE ${_link_options})
    endif()

    if(NOT "${${_feature}_TARGETS}${${_feature}_TARGETS_${_build}}" STREQUAL "")
      set(_targets ${${_feature}_TARGETS}${${_feature}_TARGETS_${_build}})
      message(STATUS "    IMPORTED TARGETS:    ${_targets}")
      copy_target_properties(${_interface_target} ${${_feature}_TARGETS} ${${_feature}_TARGETS_${_build}})
    endif()

    if(${_feature}_SPLIT_CONFIGURATION)
      # Future FIXME: change to block(PARENT_SCOPE)/endblock() (CMake 3.25)
      list(APPEND DEAL_II_TARGETS_${_build} ${_interface_target})
      set(DEAL_II_TARGETS_${_build} ${DEAL_II_TARGETS_${_build}} PARENT_SCOPE)
    else()
      # Future FIXME: change to block(PARENT_SCOPE)/endblock() (CMake 3.25)
      list(APPEND DEAL_II_TARGETS ${_interface_target})
      set(DEAL_II_TARGETS ${DEAL_II_TARGETS} PARENT_SCOPE)
    endif()

  endforeach()
endfunction()
