IF(DEAL_II_INSTALL_COMPAT_FILES)

  #
  # Build and install the old dependency_resolution executable:
  #

  ADD_EXECUTABLE(make_dependencies make_dependencies.cc)
  ADD_EXECUTABLE(report_features report_features.cc)

  ADD_DEPENDENCIES(report_features deal_ii_target_dependencies)

  INSTALL(TARGETS make_dependencies report_features
    RUNTIME
    DESTINATION common/scripts
    COMPONENT compat_files
    )

ENDIF()
