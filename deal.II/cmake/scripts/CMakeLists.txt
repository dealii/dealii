#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

ADD_EXECUTABLE(expand_instantiations expand_instantiations.cc)
IF(NOT CMAKE_CROSSCOMPILING)
  #
  # Export expand_instantiations for use in cross compilation:
  #
  EXPORT(TARGETS expand_instantiations
    FILE ${CMAKE_CURRENT_BINARY_DIR}/importExecutables.cmake
    NAMESPACE native- # Prefix a namespace
    )
ENDIF()

IF(DEAL_II_COMPONENT_COMPAT_FILES)
  #
  # Build and install the old dependency_resolution and report_features
  # executables:
  #
  INCLUDE_DIRECTORIES(
    BEFORE # Ensure deal.II include directories come first
    ${CMAKE_BINARY_DIR}/include/
    ${CMAKE_SOURCE_DIR}/include/
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
  TARGET_LINK_LIBRARIES(report_features ${DEAL_II_EXTERNAL_LIBRARIES})
  SET_TARGET_PROPERTIES(report_features
    PROPERTIES
    LINK_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}"
    COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS}"
    )

  IF(NOT CMAKE_CROSSCOMPILING)
    ADD_CUSTOM_TARGET(run_report_features COMMAND report_features)
  ENDIF()

  INSTALL(TARGETS make_dependencies report_features expand_instantiations
    EXPORT ${DEAL_II_PROJECT_CONFIG_NAME}Targets
    DESTINATION ${DEAL_II_COMMON_RELDIR}/scripts
    COMPONENT compat_files
    )
ENDIF()
