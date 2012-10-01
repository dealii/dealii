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
# Try to find SC
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(P4EST_DIR "$ENV{P4EST_DIR}")
SET_IF_EMPTY(SC_DIR "$ENV{SC_DIR}")

FIND_PATH(SC_INCLUDE_DIR sc.h
  HINTS
    ${SC_DIR}
    ${P4EST_DIR}
  PATH_SUFFIXES
    sc include/p4est include src sc/src
  )


FIND_LIBRARY(SC_LIBRARY
  NAMES sc
  HINTS
    ${SC_DIR}
    ${P4EST_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib src sc/src
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SC DEFAULT_MSG SC_LIBRARY SC_INCLUDE_DIR)

IF(SC_FOUND)
  MARK_AS_ADVANCED(
    SC_LIBRARY
    SC_INCLUDE_DIR
    SC_DIR
  )
ELSE()
  SET(SC_DIR "" CACHE STRING
    "An optional hint to an sc installation/directory"
    )
ENDIF()

