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

INCLUDE(FindPackageHandleStandardArgs)

#
# ARPACK needs LAPACK and BLAS as dependency, search for them with the help
# of the LAPACK find module:
#
# TODO: ARPACK and mpi...
#
FIND_PACKAGE(LAPACK)

FIND_LIBRARY(ARPACK_LIBRARY
  NAMES arpack
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

SET(ARPACK_LIBRARIES
  ${ARPACK_LIBRARY}
  ${LAPACK_LIBRARIES}
  )

SET(ARPACK_LINKER_FLAGS
  ${LAPACK_LINKER_FLAGS}
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(ARPACK DEFAULT_MSG
  ARPACK_LIBRARY
  LAPACK_FOUND
  )

IF(ARPACK_FOUND)
  MARK_AS_ADVANCED(
    lapack_LIBRARY
    atlas_LIBRARY
    blas_LIBRARY
    ARPACK_LIBRARY
  )
ENDIF()

