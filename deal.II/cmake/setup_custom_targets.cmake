## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
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
# Add convenience targets that build and install only a specific component:
#
#   library
#   compat_files
#   documentation
#   examples
#   mesh_converter
#   parameter_gui
#

# The library can always be installed ;-)
ADD_CUSTOM_TARGET(library
  COMMAND ${CMAKE_COMMAND}
    -DCOMPONENT="library" -P cmake_install.cmake
  COMMENT "Build and install component \"library\"."
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  )

FOREACH(_component compat_files documentation examples mesh_converter parameter_gui)
  STRING(TOUPPER "${_component}" _component_uppercase)
  IF(DEAL_II_COMPONENT_${_component_uppercase})
    ADD_CUSTOM_TARGET(${_component}
      COMMAND ${CMAKE_COMMAND}
        -DCOMPONENT="${_component}" -P cmake_install.cmake
      COMMENT "Build and install component \"${_component}\"."
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
  ELSE()
    STRING(TOUPPER ${_component} _componentuppercase)
    ADD_CUSTOM_TARGET(${_component}
      COMMAND
           ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo '***************************************************************************'
        && ${CMAKE_COMMAND} -E echo "**  Error: Could not build and install disabled component \"${_component}\"."
        && ${CMAKE_COMMAND} -E echo "**  Please reconfigure with -DDEAL_II_COMPONENT_${_componentuppercase}=ON"
        && ${CMAKE_COMMAND} -E echo '***************************************************************************'
        && ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo ''
        && false
      )
  ENDIF()
ENDFOREACH()

#
# Provide an "info" target to print a help message:
#

FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
"MESSAGE(
\"###
#
#  The following targets are available (invoke by $ make <target>):
#
#    all            - compiles the library and all enabled components
#    clean          - removes all generated files
#    install        - installs into CMAKE_INSTALL_PREFIX
#    help           - prints a list of valid top level targets
#    info           - prints this help message
#
#    edit_cache     - runs ccmake for changing (cached) configuration variables
#                     and reruns the configure and generate phases of CMake
#    rebuild_cache  - reruns the configure and generate phases of CMake
#
#    compat_files   - builds and installs the 'compat_files' component
#    documentation  - builds and installs the 'documentation' component
#    examples       - builds and installs the 'examples' component
#    library        - builds and installs the 'library' component
#    mesh_converter - builds and installs the 'mesh_converter' component
#    parameter_gui  - builds and installs the 'parameter_gui' component
#
#    test           - runs a minimal set of tests
#
#    setup_tests    - sets up the testsuite subprojects
#    clean_tests    - runs the 'clean' target in every testsuite subproject
#    prune_tests    - removes all testsuite subprojects
#
###\")"
  )
ADD_CUSTOM_TARGET(info
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
  )
