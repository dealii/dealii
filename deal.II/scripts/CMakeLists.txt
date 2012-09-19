IF(DEAL_II_INSTALL_COMPAT_FILES)

  #
  # Build and install the old dependency_resolution executable:
  #

  ADD_EXECUTABLE(make_dependencies make_dependencies.cc)

  INSTALL(TARGETS make_dependencies
    RUNTIME
    DESTINATION common/scripts
    COMPONENT compat_files
    )

ENDIF()
