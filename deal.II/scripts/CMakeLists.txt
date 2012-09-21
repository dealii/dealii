IF(DEAL_II_INSTALL_COMPAT_FILES)

  #
  # Build and install the old dependency_resolution executable:
  #

  INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/include/
    ${CMAKE_BINARY_DIR}/include/
    )

  ADD_EXECUTABLE(make_dependencies make_dependencies.cc)

  ADD_EXECUTABLE(report_features report_features.cc)
  TARGET_LINK_LIBRARIES(report_features ${DEAL_II_EXTERNAL_LIBRARIES})

  ADD_DEPENDENCIES(report_features deal_ii_target_dependencies)

  INSTALL(TARGETS make_dependencies report_features
    RUNTIME
    DESTINATION common/scripts
    COMPONENT compat_files
    )

ENDIF()
