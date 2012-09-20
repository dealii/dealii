
#
# Set up all necessary bits for the documentation
#

IF(DEAL_II_WITH_DOXYGEN)

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
    FILES_MATCHING REGEX "^.*(html|css|136|eps|fig|png)$"  # TODO! cmapx, dot
    )

  SET(doc_files_to_configure
    title.html.in
    doxygen/tutorial/toc.html.in
    doxygen/header.136.in
    doxygen/header.tex.in
    doxygen/options.dox.in
    doxygen/Makefile.in
    doxygen/options.136.in
    doxygen/header.html.in
    )

  FOREACH(infile ${doc_files_to_configure})
    STRING(REGEX REPLACE "\\.in$" "" outfile "${infile}" )
    CONFIGURE_FILE(
      ${CMAKE_CURRENT_SOURCE_DIR}/${infile}
      ${CMAKE_CURRENT_BINARY_DIR}/${outfile}
      )
    INSTALL(FILES
      ${CMAKE_CURRENT_BINARY_DIR}/${outfile}
      DESTINATION ${DEAL_II_DOCUMENTATION_RELDIR}
      COMPONENT documentation
      )
  ENDFOREACH()

ENDIF(DEAL_II_WITH_DOXYGEN)

