## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Configuration for the Verdict library:
#


MACRO(FEATURE_VERDICT_CONFIGURE_COMMON)

ENDMACRO()


MACRO(FEATURE_VERDICT_CONFIGURE_BUNDLED)
  SET(VERDICT_BUNDLED_INCLUDE_DIRS ${VERDICT_FOLDER})

  FEATURE_VERDICT_CONFIGURE_COMMON()
ENDMACRO()

MACRO(FEATURE_VERDICT_FIND_EXTERNAL var)
  FIND_PACKAGE(VERDICT)

  IF(VERDICT_FOUND)
    SET(${var} TRUE)
  ENDIF()
ENDMACRO()


MACRO(FEATURE_VERDICT_CONFIGURE_EXTERNAL)
  FEATURE_VERDICT_CONFIGURE_COMMON()
ENDMACRO()


CONFIGURE_FEATURE(VERDICT)
