#
# Configuration for doxygen
#

MACRO(FEATURE_DOXYGEN_FIND_EXTERNAL var)

  FIND_PACKAGE(Doxygen)

  #
  # We use doxygen and dot:
  #
  IF(DOXYGEN_FOUND AND DOXYGEN_DOT_FOUND)
    SET(${var} TRUE)
  ENDIF()

ENDMACRO()


MACRO(FEATURE_DOXYGEN_CONFIGURE_EXTERNAL var)

  #
  # The FindDoxygen defines for us:
  #
  #   DOXYGEN_EXECUTABLE     = The path to the doxygen command.
  #   DOXYGEN_FOUND          = Was Doxygen found or not?
  #   DOXYGEN_VERSION        = The version reported by doxygen --version
  #
  #   DOXYGEN_DOT_EXECUTABLE = The path to the dot program used by doxygen.
  #   DOXYGEN_DOT_FOUND      = Was Dot found or not?
  #   DOXYGEN_DOT_PATH       = The path to dot not including the executable
  #
  # Use these variables to set up a custom target:
  #

  SET(${var} TRUE)

ENDMACRO()


SET(FEATURE_DOXYGEN_CUSTOM_ERROR_MESSAGE TRUE)


MACRO(FEATURE_DOXYGEN_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "\n"
    "Could not find the doxygen package!\n\n"
    "Please ensure that doxygen and dot are installed on your computer.\n"
    "If the packages are not at a default location, either provide some hints\n"
    "for the autodetection, or set the relevant variables by hand in ccmake.\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(DOXYGEN)

