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
# Try to find the ARPACK library
#
# This module exports
#
#   ARPACK_LIBRARIES
#   ARPACK_LINKER_FLAGS
#

#
# TODO: ARPACK and mpi...
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(ARPACK_DIR "$ENV{ARPACK_DIR}")

#
# ARPACK needs LAPACK and BLAS as dependencies:
#
FIND_PACKAGE(LAPACK)

FIND_LIBRARY(ARPACK_LIBRARY
  NAMES arpack
  HINTS
    ${ARPACK_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(ARPACK DEFAULT_MSG
  ARPACK_LIBRARY
  LAPACK_FOUND
  )

IF(ARPACK_FOUND)
  SET(ARPACK_LIBRARIES
    ${ARPACK_LIBRARY}
    ${LAPACK_LIBRARIES}
    )
  SET(ARPACK_LINKER_FLAGS
    ${LAPACK_LINKER_FLAGS}
    )

  MARK_AS_ADVANCED(
    lapack_LIBRARY
    atlas_LIBRARY
    blas_LIBRARY
    ARPACK_LIBRARY
    ARPACK_DIR
  )
ELSE()
  SET(ARPACK_DIR "" CACHE PATH
    "An optional hint to an ARPACK installation"
    )
ENDIF()

