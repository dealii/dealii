## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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

ADD_EXECUTABLE(expand_instantiations_exe expand_instantiations.cc)
SET_TARGET_PROPERTIES(expand_instantiations_exe
  PROPERTIES OUTPUT_NAME expand_instantiations
  )

IF(DEAL_II_COMPONENT_COMPAT_FILES)
  MESSAGE(STATUS "Setting up compat_files")

  #
  # Also build and install the old dependency_resolution executable:
  #

  ADD_EXECUTABLE(make_dependencies make_dependencies.cc)

  ADD_DEPENDENCIES(compat_files
    expand_instantiations_exe
    make_dependencies
    )

  SET_TARGET_PROPERTIES(expand_instantiations_exe make_dependencies
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${DEAL_II_COMMON_RELDIR}/scripts
    )
  INSTALL(TARGETS expand_instantiations_exe make_dependencies
    DESTINATION ${DEAL_II_COMMON_RELDIR}/scripts
    COMPONENT compat_files
    )
  MESSAGE(STATUS "Setting up compat_files - Done")
ENDIF()
