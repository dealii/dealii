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
# Configuration for the netcdf library:
#


MACRO(FEATURE_METIS_FIND_EXTERNAL var)
  FIND_PACKAGE(METIS)

  IF(METIS_FOUND AND METIS_MAJOR GREATER 4)
    SET(${var} TRUE)
  ELSE()
    MESSAGE(WARNING "\n"
      "Could not find a sufficient modern metis installation: "
      "Version 5.x required!\n\n"
      )
  ENDIF()
ENDMACRO()


MACRO(FEATURE_METIS_CONFIGURE_EXTERNAL var)
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${METIS_LIBRARY})

  SET(DEAL_II_USE_METIS TRUE)

  SET(${var} TRUE)
ENDMACRO()


SET(FEATURE_METIS_CUSTOM_ERROR_MESSAGE TRUE)


MACRO(FEATURE_METIS_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find the metis library!\n\n"
    "Please ensure that the metis library version 5.0 or newer is installed on your computer.\n"
    "If the library is not at a default location, either provide some hints\n"
    "for the autodetection:\n"
    "    $ METIS_DIR=\"...\" cmake <...>\n"
    "    $ ccmake -DMETIS_DIR=\"...\" cmake <...>\n"
    "or set the relevant variables by hand in ccmake.\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(METIS)

