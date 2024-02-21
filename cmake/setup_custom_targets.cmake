## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2022 by the deal.II authors
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
# Add convenience targets that build and install only a specific component:
#
#   library
#   documentation
#   examples
#


if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local")
  #
  # In case that CMAKE_INSTALL_PREFIX wasn't set, we assume that the user
  # doesn't actually want to install but just use deal.II in the build
  # directory. In this case, do not add the "install" phase to the
  # convenience targets.
  #
  macro(_add_custom_target _name)
    add_custom_target(${_name})
  endmacro()

  # Print precise information about the convenience targets:
  set(_description_string "build")
else()
  macro(_add_custom_target _name)
    add_custom_target(${_name}
      COMMAND ${CMAKE_COMMAND}
        -DCOMPONENT="${_name}" -P cmake_install.cmake
      COMMENT "Build and install component \"library\"."
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
  endmacro()

  # Print precise information about the convenience targets:
  set(_description_string "build and install")
endif()

# The library can always be compiled and/or installed unconditionally ;-)
_add_custom_target(library)

foreach(_component documentation examples python_bindings)
  string(TOUPPER "${_component}" _component_uppercase)
  if(DEAL_II_COMPONENT_${_component_uppercase})
    _add_custom_target(${_component})
  else()
    string(TOUPPER ${_component} _componentuppercase)

    set(_error_description_message
      "Error: Could not ${_description_string} disabled component ${_component}.")
    decorate_with_stars(${_error_description_message}
      _decorated_error_description_message)

    set(_reconfiguration_help_message
      "Please reconfigure with -DDEAL_II_COMPONENT_${_componentuppercase}=ON")
    decorate_with_stars(${_reconfiguration_help_message}
      _decorated_reconfiguration_help_message)

    add_custom_target(${_component}
      COMMAND
           ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo '***************************************************************************'
        && ${CMAKE_COMMAND} -E echo "${_decorated_error_description_message}"
        && ${CMAKE_COMMAND} -E echo "${_decorated_reconfiguration_help_message}"
        && ${CMAKE_COMMAND} -E echo '***************************************************************************'
        && ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo ''
        && false
      )
  endif()
endforeach()

if(NOT DEAL_II_COMPONENT_PACKAGE)
  add_custom_target(package
    COMMAND
         ${CMAKE_COMMAND} -E echo ''
      && ${CMAKE_COMMAND} -E echo ''
      && ${CMAKE_COMMAND} -E echo '***************************************************************************'
      && ${CMAKE_COMMAND} -E echo '**  Error: Could not generate binary package. The component is disabled. **'
      && ${CMAKE_COMMAND} -E echo '**        Please reconfigure with -DDEAL_II_COMPONENT_PACKAGE=ON         **'
      && ${CMAKE_COMMAND} -E echo '***************************************************************************'
      && ${CMAKE_COMMAND} -E echo ''
      && ${CMAKE_COMMAND} -E echo ''
      && false
    )
endif()

#
# Provide an indentation target for indenting uncommitted changes and changes on
# the current feature branch
#
add_custom_target(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./contrib/utilities/indent
  COMMENT "Indenting recently changed files in the deal.II directories"
  )

#
# Provide "indent" target for indenting all headers and source files
#
add_custom_target(indent-all
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./contrib/utilities/indent-all
  COMMENT "Indenting all files in the deal.II directories"
  )


#
# Provide an "info" target to print a help message:
#
if(CMAKE_GENERATOR MATCHES "Ninja")
  set(_make_command "ninja")
else()
  set(_make_command "make")
endif()

file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
"message(
\"###
#
#  The following targets are available (invoke by $ ${_make_command} <target>):
#
#    all            - compile the library and all enabled components
#    clean          - remove all generated files
#    install        - install into CMAKE_INSTALL_PREFIX
#
#    info           - print this help message
#    help           - print a list of valid top level targets
#
#    edit_cache     - run ccmake for changing (cached) configuration variables
#                     and reruns the configure and generate phases of CMake
#    rebuild_cache  - rerun the configure and generate phases of CMake
#
#    documentation  - ${_description_string} component 'documentation'
#    examples       - ${_description_string} component 'examples'
#    library        - ${_description_string} component 'library'
#    package        - build binary package
#
#    test           - run a minimal set of tests
#
#    setup_tests    - set up testsuite subprojects
#    prune_tests    - remove all testsuite subprojects
#
#    indent         - indent all headers and source files that changed since the
#                     last commit to master, including untracked ones
#    indent-all     - indent all headers and source files
")

#
# Provide "relocate" target to run install_name_tool on all external libraries
# under ${DEAL_II_CPACK_EXTERNAL_LIBS_TREE}
#
if(CMAKE_SYSTEM_NAME MATCHES "Darwin" AND
  NOT "${DEAL_II_CPACK_EXTERNAL_LIBS_TREE}" STREQUAL "")
  add_custom_target(relocate
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMAND ./contrib/utilities/relocate_libraries.py
    COMMENT "Running install_name_tool under ${DEAL_II_CPACK_EXTERNAL_LIBS_TREE}"
    )
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
  "#
#    relocate       - fix RPATH for external libraries, if packaging was requested
"
   )
endif()

file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
"#
###\")"
)

add_custom_target(info
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
  )


#
# Provide a target to build all .inst files
#
add_custom_target(expand_all_instantiations)
