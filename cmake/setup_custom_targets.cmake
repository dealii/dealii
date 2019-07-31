## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2018 by the deal.II authors
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
# Add convenience targets that build and install only a specific component:
#
#   library
#   documentation
#   examples
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

  # Print precise information about the convenience targets:
  SET(_description_string "build")
ELSE()
  MACRO(_add_custom_target _name)
    ADD_CUSTOM_TARGET(${_name}
      COMMAND ${CMAKE_COMMAND}
        -DCOMPONENT="${_name}" -P cmake_install.cmake
      COMMENT "Build and install component \"library\"."
      WORKING_DIRECTORY ${deal.II_BINARY_DIR}
      )
  ENDMACRO()

  # Print precise information about the convenience targets:
  SET(_description_string "build and install")
ENDIF()

# The library can always be compiled and/or installed unconditionally ;-)
_add_custom_target(library)

FOREACH(_component documentation examples python_bindings)
  STRING(TOUPPER "${_component}" _component_uppercase)
  IF(DEAL_II_COMPONENT_${_component_uppercase})
    _add_custom_target(${_component})
  ELSE()
    STRING(TOUPPER ${_component} _componentuppercase)

    SET(_error_description_message
      "Error: Could not ${_description_string} disabled component ${_component}.")
    DECORATE_WITH_STARS(${_error_description_message}
      _decorated_error_description_message)

    SET(_reconfiguration_help_message
      "Please reconfigure with -DDEAL_II_COMPONENT_${_componentuppercase}=ON")
    DECORATE_WITH_STARS(${_reconfiguration_help_message}
      _decorated_reconfiguration_help_message)

    ADD_CUSTOM_TARGET(${_component}
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
  ENDIF()
ENDFOREACH()

IF(NOT DEAL_II_COMPONENT_PACKAGE)
  ADD_CUSTOM_TARGET(package
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
ENDIF()

#
# Provide an indentation target for indenting uncommitted changes and changes on
# the current feature branch
#
ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${deal.II_SOURCE_DIR}
  COMMAND ./contrib/utilities/indent
  COMMENT "Indenting recently changed files in the deal.II directories"
  )

#
# Provide "indent" target for indenting all headers and source files
#
ADD_CUSTOM_TARGET(indent-all
  WORKING_DIRECTORY ${deal.II_SOURCE_DIR}
  COMMAND ./contrib/utilities/indent-all
  COMMENT "Indenting all files in the deal.II directories"
  )


#
# Provide an "info" target to print a help message:
#
IF(CMAKE_GENERATOR MATCHES "Ninja")
  SET(_make_command "ninja")
ELSE()
  SET(_make_command "make")
ENDIF()

FILE(WRITE ${deal.II_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
"MESSAGE(
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
IF(CMAKE_SYSTEM_NAME MATCHES "Darwin" AND
  NOT "${DEAL_II_CPACK_EXTERNAL_LIBS_TREE}" STREQUAL "")
  ADD_CUSTOM_TARGET(relocate
    WORKING_DIRECTORY ${deal.II_SOURCE_DIR}
    COMMAND ./contrib/utilities/relocate_libraries.py
    COMMENT "Running install_name_tool under ${DEAL_II_CPACK_EXTERNAL_LIBS_TREE}"
    )
  FILE(APPEND ${deal.II_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
  "#
#    relocate       - fix RPATH for external libraries, if packaging was requested
"
   )
ENDIF()

FILE(APPEND ${deal.II_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
"#
###\")"
)

ADD_CUSTOM_TARGET(info
  COMMAND ${CMAKE_COMMAND} -P ${deal.II_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/print_info.cmake
  )


#
# Provide a target to build all .inst files
#
ADD_CUSTOM_TARGET(expand_all_instantiations)
