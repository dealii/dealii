## ---------------------------------------------------------------------
## $Id$
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

IF(NOT CMAKE_CROSSCOMPILING)
  ADD_EXECUTABLE(expand_instantiations expand_instantiations.cc)

  IF(DEAL_II_COMPONENT_COMPAT_FILES)
    SET_TARGET_PROPERTIES(expand_instantiations
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${DEAL_II_COMMON_RELDIR}/scripts
      )
    INSTALL(TARGETS expand_instantiations
      DESTINATION ${DEAL_II_COMMON_RELDIR}/scripts
      COMPONENT compat_files
      )
  ENDIF()

  EXPORT(TARGETS expand_instantiations
    FILE ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Executables.cmake
    )
ENDIF()

IF(DEAL_II_COMPONENT_COMPAT_FILES)
  #
  # Build and install the old dependency_resolution and report_features
  # executables:
  #

  INCLUDE_DIRECTORIES(
    BEFORE # Ensure deal.II include directories come first
    #
    # Reverse order due to BEFORE:
    #
    ${CMAKE_SOURCE_DIR}/include/
    ${CMAKE_BINARY_DIR}/include/
    )

  ADD_EXECUTABLE(make_dependencies make_dependencies.cc)

  IF(BUILD_SHARED_LIBS)
    SET(_shared_libs "true")
  ELSE()
    SET(_shared_libs "false")
  ENDIF()
  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/report_features.cc.in
    ${CMAKE_CURRENT_BINARY_DIR}/report_features.cc
    )

  ADD_EXECUTABLE(report_features
    ${CMAKE_CURRENT_BINARY_DIR}/report_features.cc
    )

  #
  # Link report_features in the same way as the library/final executable:
  #
  TARGET_LINK_LIBRARIES(report_features ${DEAL_II_EXTERNAL_LIBRARIES})
  SET_TARGET_PROPERTIES(report_features
    PROPERTIES
    LINK_FLAGS "${DEAL_II_LINKER_FLAGS}"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS}"
    )
  IF(DEAL_II_STATIC_EXECUTABLE)
    SET_PROPERTY(TARGET report_features PROPERTY
      LINK_SEARCH_END_STATIC TRUE
      )
  ENDIF()

  ADD_DEPENDENCIES(compat_files
    expand_instantiations
    make_dependencies
    report_features
    )

  IF(NOT CMAKE_CROSSCOMPILING)
    ADD_CUSTOM_TARGET(run_report_features COMMAND report_features)
  ENDIF()

  SET_TARGET_PROPERTIES(make_dependencies report_features
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${DEAL_II_COMMON_RELDIR}/scripts
    )

  INSTALL(TARGETS make_dependencies report_features
    DESTINATION ${DEAL_II_COMMON_RELDIR}/scripts
    COMPONENT compat_files
    )
  EXPORT(TARGETS make_dependencies report_features
    FILE ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Executables.cmake
    APPEND
    )
ENDIF()
