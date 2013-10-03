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

FOREACH(_component library)
  ADD_CUSTOM_TARGET(${_component}
    COMMAND ${CMAKE_COMMAND}
      -DCOMPONENT="${_component}" -P cmake_install.cmake
    COMMENT "Build and install component \"${_component}\"."
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
ENDFOREACH()

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

    ADD_CUSTOM_TARGET(${_component}
      COMMAND ${CMAKE_COMMAND} -E echo "Error: Could not build and install disabled component \"${_component}\"."
        && ${CMAKE_COMMAND} -E echo "Please reconfigure with -DDEAL_II_COMPONENT_${_component}=yes"
        && false
      )

  ENDIF()
ENDFOREACH()

