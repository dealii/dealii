#
# Configuration for the umfpack and amd libraries:
#


SET(FEATURE_UMFPACK_DEPENDS
  DEAL_II_WITH_BLAS
  DEAL_II_WITH_LAPACK
  )


MACRO(FEATURE_UMFPACK_FIND_EXTERNAL var)
  FIND_PACKAGE(UMFPACK)
  FIND_PACKAGE(AMD)

  IF(UMFPACK_FOUND AND AMD_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_UMFPACK_CONFIGURE_EXTERNAL var)
  INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES
    ${UMFPACK_LIBRARY} ${AMD_LIBRARY}
    )

  SET(HAVE_LIBUMFPACK TRUE)

  SET(${var} TRUE)
ENDMACRO()


SET(FEATURE_UMFPACK_HAVE_CONTRIB TRUE)


MACRO(FEATURE_UMFPACK_CONFIGURE_CONTRIB var)
  #
  # Add umfpack and amd directly to the object files of deal.II
  #

  SET(umfpack_folder "${CMAKE_SOURCE_DIR}/contrib/umfpack")

  INCLUDE_DIRECTORIES(
    ${umfpack_folder}/UMFPACK/Include
    ${umfpack_folder}/AMD/Include
    )

  ADD_SUBDIRECTORY(${umfpack_folder}/UMFPACK/Source)
  ADD_SUBDIRECTORY(${umfpack_folder}/AMD/Source)

  SET(HAVE_LIBUMFPACK TRUE)

  SET(${var} TRUE)
ENDMACRO()


SET(FEATURE_UMFPACK_CUSTOM_ERROR_MESSAGE TRUE)


MACRO(FEATURE_UMFPACK_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "\n"
    "Could not find the umfpack and amd libraries!\n"
    "Please ensure that the libraries are installed on your computer.\n"
    "If the libraries are not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "    $ UMFPACK_DIR=\"...\" cmake <...>\n"
    "    $ ccmake -DUMFPACK_DIR=\"...\" cmake <...>\n"
    "or set the relevant variables by hand in ccmake.\n"
    "Alternatively you may choose to compile the bundled contrib libraries\n"
    "by setting DEAL_II_ALLOW_CONTRIB=on or DEAL_II_FORCE_CONTRIB_UMFPACK=on.\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(UMFPACK)

