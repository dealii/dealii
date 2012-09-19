
ADD_EXECUTABLE(expand_instantiations expand_instantiations.cc)

IF(DEAL_II_INSTALL_COMPAT_FILES)
  #
  # Install the expand_instantiations to its old place under common/scripts
  #
  INSTALL(TARGETS expand_instantiations
    RUNTIME
    DESTINATION common/scripts
    COMPONENT compat_files
    )

ENDIF()
