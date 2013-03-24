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
# Try to find the Threading Building Blocks library
#
# This module exports
#
#   TBB_LIBRARIES
#   TBB_INCLUDE_DIRS
#   TBB_WITH_DEBUGLIB
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(TBB_DIR "$ENV{TBB_DIR}")

FIND_PATH(TBB_INCLUDE_DIR tbb/parallel_reduce.h
  HINTS
    ${TBB_DIR}
  PATH_SUFFIXES include include/tbb tbb)

FIND_LIBRARY(TBB_LIBRARY
  NAMES tbb
  HINTS
    ${TBB_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

FIND_LIBRARY(TBB_DEBUG_LIBRARY
  NAMES tbb_debug
  HINTS
    ${TBB_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(TBB DEFAULT_MSG
  TBB_LIBRARY
  TBB_INCLUDE_DIR
  )

IF(TBB_FOUND)
  MARK_AS_ADVANCED(TBB_LIBRARY TBB_DEBUG_LIBRARY TBB_INCLUDE_DIR)

  IF(NOT TBB_DEBUG_LIBRARY MATCHES "-NOTFOUND")
    SET(TBB_WITH_DEBUGLIB TRUE)
    SET(TBB_LIBRARIES debug ${TBB_LIBARY} optimized ${TBB_LIBRARY})
  ELSE()
    SET(TBB_LIBRARIES ${TBB_LIBARY})
  ENDIF()

  SET(TBB_INCLUDE_DIRS ${TBB_INCLUDE_DIR})

ELSE()

  SET(TBB_DIR "" CACHE PATH
    "An optional hint to a TBB installation"
    )
ENDIF()

