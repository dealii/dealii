#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Configuration for doxygen
#

OPTION(DEAL_II_WITH_DOXYGEN
  "Build deal.II with support for doxygen and dot."
  OFF)


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

#
# Additional dependency check:
# DEAL_II_COMPONENT_DOCUMENTATION needs DEAL_II_WITH_DOXYGEN
#
IF(DEAL_II_COMPONENT_DOCUMENTATION AND NOT DEAL_II_WITH_DOXYGEN)
  IF(DEAL_II_FEATURE_AUTODETECTION)
    FEATURE_DOXYGEN_ERROR_MESSAGE()
  ELSE()
    MESSAGE(SEND_ERROR "\n"
      "DEAL_II_COMPONENT_DOCUMENTATION has unmet configuration requirements: "
      "DEAL_II_WITH_DOXYGEN required, but set to OFF!\n\n"
      )
  ENDIF()
ENDIF()

