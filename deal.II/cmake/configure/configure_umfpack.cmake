#
# Configuration for the umfpack and amd libraries:
#


SET(FEATURE_UMFPACK_DEPENDS
  # Currently, with enabled umfpack support, we also need to setup
  # LAPACK support in deal.II:
  DEAL_II_WITH_LAPACK
  )


MACRO(FEATURE_UMFPACK_FIND_EXTERNAL var)
  FIND_PACKAGE(UMFPACK)

  IF(UMFPACK_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_UMFPACK_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES(${UMFPACK_INCLUDE_DIRS})
  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${UMFPACK_LIBRARIES})
  ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${UMFPACK_LINKER_FLAGS}")

  SET(HAVE_LIBUMFPACK TRUE)

  SET(${var} TRUE)
ENDMACRO()


SET(FEATURE_UMFPACK_HAVE_BUNDLED TRUE)


MACRO(FEATURE_UMFPACK_CONFIGURE_BUNDLED var)
  #
  # DEAL_II_WITH_LAPACK will pull in an external BLAS library. So no need
  # to setup something more than bundled UMFPACK here.
  #

  SET(umfpack_folder "${CMAKE_SOURCE_DIR}/bundled/umfpack")
  INCLUDE_DIRECTORIES(
    ${umfpack_folder}/UMFPACK/Include
    ${umfpack_folder}/AMD/Include
    )

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
    "Alternatively you may choose to compile the bundled libraries\n"
    "by setting DEAL_II_ALLOW_BUNDLED=ON or DEAL_II_FORCE_BUNDLED_UMFPACK=ON.\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(UMFPACK)

