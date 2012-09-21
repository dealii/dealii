CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
  )


IF(DEAL_II_INSTALL_PROJECT_CONFIG)
  #
  # Configure and install the cmake project configuration:
  #

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/deal.IIConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/deal.IIConfig.cmake
    )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/deal.IIConfigVersion.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/deal.IIConfigVersion.cmake
    @ONLY
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/deal.IIConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/deal.IIConfigVersion.cmake
    DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
    COMPONENT project_config
    )
ENDIF()


IF(DEAL_II_INSTALL_COMPAT_FILES)
  #
  # Configure and install the old Makefile configuration:
  #

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.in
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Version.in
    ${CMAKE_CURRENT_BINARY_DIR}/Version
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    DESTINATION common
    COMPONENT compat_files
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/Version
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    COMPONENT compat_files
    )

ENDIF()
