
#
# Set up all necessary bits for the documentation
#

IF(DEAL_II_WITH_DOXYGEN)

  #
  # Generate the documentation via doxygen:
  #
  ADD_SUBDIRECTORY(doxygen)


  #
  # Install the homepage:
  #

  IF(DEAL_II_INSTALL_COMPAT_FILES)
    #
    # The good, old directory where the documentation gets installed to:
    #
    SET(DEAL_II_DOCUMENTATION_RELDIR "doc")
  ELSE()
    #
    # IF DEAL_II_INSTALL_COMPAT_FILES is not set, we assume that we have to
    # obey the file system standard, so install to usr/share/doc/deal.II.
    #
    SET(DEAL_II_DOCUMENTATION_RELDIR "usr/share/doc/deal.II/html")
  ENDIF()

  INSTALL(DIRECTORY .
    DESTINATION ${DEAL_II_DOCUMENTATION_RELDIR}
    COMPONENT documentation
    FILES_MATCHING REGEX "^.*(html|css|eps|fig|png)$"  # TODO! cmapx, dot
    )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/title.html.in
    ${CMAKE_CURRENT_BINARY_DIR}/title.html
    )
  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/title.html
    DESTINATION ${DEAL_II_DOCUMENTATION_RELDIR}
    COMPONENT documentation
    )

ENDIF(DEAL_II_WITH_DOXYGEN)

