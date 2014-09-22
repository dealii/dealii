## ---------------------------------------------------------------------
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


IF("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local")
  #
  # In case that CMAKE_INSTALL_PREFIX wasn't set, we assume that the user
  # doesn't actually want to install but just use deal.II in the build
  # directory. In this case, do not add the "install" phase to the
  # convenience targets.
  #
  MACRO(_add_custom_target _name)
    ADD_CUSTOM_TARGET(${_name})
  ENDMACRO()

  # Print precise informations about the convenience targets:
  SET(_description_string "build")
ELSE()
  MACRO(_add_custom_target _name)
    ADD_CUSTOM_TARGET(${_name}
      COMMAND ${CMAKE_COMMAND}
        -DCOMPONENT="${_name}" -P cmake_install.cmake
      COMMENT "Build and install component \"library\"."
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      )
  ENDMACRO()

  # Print precise informations about the convenience targets:
  SET(_description_string "build and install")
ENDIF()

# The library can always be compiled and/or installed unconditionally ;-)
_add_custom_target(library)

FOREACH(_component compat_files documentation examples mesh_converter parameter_gui)
  STRING(TOUPPER "${_component}" _component_uppercase)
  IF(DEAL_II_COMPONENT_${_component_uppercase})
    _add_custom_target(${_component})
  ELSE()
    STRING(TOUPPER ${_component} _componentuppercase)
    ADD_CUSTOM_TARGET(${_component}
      COMMAND
           ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo ''
        && ${CMAKE_COMMAND} -E echo '***************************************************************************'
        && ${CMAKE_COMMAND} -E echo "**  Error: Could not ${_description_string} disabled component \"${_component}\"."
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
IF(CMAKE_GENERATOR MATCHES "Ninja")
  SET(_make_command "ninja")
ELSE()
  SET(_make_command "make")
ENDIF()

FILE(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
"MESSAGE(
\"###
#
#  The following targets are available (invoke by $ ${_make_command} <target>):
#
#    all            - compile the library and all enabled components
#    clean          - remove all generated files
#    install        - install into CMAKE_INSTALL_PREFIX
#    info           - print this help message
#    help           - print a list of valid top level targets
#
#    edit_cache     - run ccmake for changing (cached) configuration variables
#                     and reruns the configure and generate phases of CMake
#    rebuild_cache  - rerun the configure and generate phases of CMake
#
#    compat_files   - ${_description_string} component 'compat_files'
#    documentation  - ${_description_string} component 'documentation'
#    examples       - ${_description_string} component 'examples'
#    library        - ${_description_string} component 'library'
#    mesh_converter - ${_description_string} component 'mesh_converter'
#    parameter_gui  - ${_description_string} component 'parameter_gui'
#
#    test           - run a minimal set of tests
#
#    setup_tests    - set up testsuite subprojects
#    regen_tests    - rerun configure stage in every testsuite subproject
#    clean_tests    - run the 'clean' target in every testsuite subproject
#    prune_tests    - remove all testsuite subprojects
#
###\")"
  )

ADD_CUSTOM_TARGET(info
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
  )
