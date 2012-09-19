CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
  )


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
