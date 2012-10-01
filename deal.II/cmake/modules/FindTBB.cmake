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

INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(TBB_INCLUDE_DIR tbb/parallel_reduce.h
)

FIND_LIBRARY(TBB_LIBRARY
  NAMES tbb
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)

FIND_LIBRARY(TBB_DEBUG_LIBRARY
  NAMES tbb_debug
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(TBB DEFAULT_MSG TBB_LIBRARY TBB_INCLUDE_DIR)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TBB_DEBUG DEFAULT_MSG TBB_DEBUG_LIBRARY TBB_INCLUDE_DIR)

IF(TBB_FOUND)
  MARK_AS_ADVANCED(
    TBB_LIBRARY
    TBB_DEBUG_LIBRARY
    TBB_INCLUDE_DIR
  )
ENDIF()

