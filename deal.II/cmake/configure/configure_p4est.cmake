#
# Configuration for the p4est and sc libraries:
#

OPTION(DEAL_II_WITH_P4EST
  "Build deal.II with support for p4est."
  OFF)


SET(FEATURE_P4EST_DEPENDS DEAL_II_WITH_MPI)


MACRO(FEATURE_P4EST_FIND_EXTERNAL var)

  FIND_PACKAGE(P4EST)

  IF(P4EST_FOUND)
    #
    # Check whether p4est supports mpi:
    #
    IF(NOT P4EST_WITH_MPI)
      MESSAGE(WARNING "\n"
        "Could not find a sufficient p4est installation: "
        "P4est has to be configured with MPI support enabled.\n\n"
        )
    ELSE()
      SET(${var} TRUE)
    ENDIF()
  ENDIF()

ENDMACRO()


MACRO(FEATURE_P4EST_CONFIGURE_EXTERNAL var)

  INCLUDE_DIRECTORIES(${P4EST_INCLUDE_DIRS})

  # The user has to know the location of the p4est headers as well:
  LIST(APPEND DEAL_II_USER_INCLUDE_DIRS ${P4EST_INCLUDE_DIRS})

  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${P4EST_LIBRARIES})

  SET(DEAL_II_USE_P4EST TRUE)

  SET(${var} TRUE)
ENDMACRO()


SET(FEATURE_P4EST_CUSTOM_ERROR_MESSAGE TRUE)


MACRO(FEATURE_P4EST_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "\n"
    "Could not find the p4est and sc libraries!\n\n"
    "Please ensure that the libraries are installed on your computer.\n"
    "If the libraries are not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "    $ P4EST_DIR=\"...\" cmake <...>\n"
    "    $ ccmake -DP4EST_DIR=\"...\" cmake <...>\n"
    "or set the relevant variables by hand in ccmake.\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(P4EST)

