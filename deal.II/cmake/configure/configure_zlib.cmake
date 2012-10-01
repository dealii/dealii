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
# Configuration for the zlib library:
#

OPTION(DEAL_II_WITH_ZLIB
  "Build deal.II with support for zlib."
  OFF)

MACRO(FEATURE_ZLIB_FIND_EXTERNAL var)
  FIND_PACKAGE(ZLIB)

  IF(ZLIB_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_ZLIB_CONFIGURE_EXTERNAL var)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${ZLIB_LIBRARIES})
  SET(HAVE_LIBZ TRUE)

  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(ZLIB)

