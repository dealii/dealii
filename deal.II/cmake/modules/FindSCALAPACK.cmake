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
# Try to find the SCALAPACK library
#
# Used as a helper module for FindMUMPS.cmake
#
# This module exports
#
#   SCALAPACK_LIBRARIES
#   SCALAPACK_LINKER_FLAGS
#

SET_IF_EMPTY(SCALAPACK_DIR "$ENV{SCALAPACK_DIR}")
SET_IF_EMPTY(BLACS_DIR "$ENV{BLACS_DIR}")

INCLUDE(FindPackageHandleStandardArgs)

FIND_LIBRARY(SCALAPACK_LIBRARY NAMES scalapack
  HINTS
    ${SCALAPACK_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
  )

#
# SCALAPACK needs LAPACK and BLAS as dependency, search for them with the help
# of the LAPACK find module:
#
FIND_PACKAGE(LAPACK)

#
# Well, depending on the version of scalapack and the distribution it might
# be necessary to search for blacs, too. So we do this in a very
# probabilistic way...
#
FOREACH(_lib blacs blacsCinit blacsF77init)
  STRING(TOUPPER "${_lib}" _lib_upper)
  FIND_LIBRARY(${_lib_upper}_LIBRARY
    NAMES ${_lib} ${_lib}_MPI-LINUX-0
    HINTS
      ${BLACS_DIR}
      ${SCALAPACK_DIR}
      ${SCALAPACK_DIR}/../blacs/
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib LIB
  )
  IF(NOT ${_lib_upper}_LIBRARY MATCHES "-NOTFOUND")
    LIST(APPEND BLACS_LIBRARIES
      ${${_lib_upper}_LIBRARY}
      )
  ENDIF()
ENDFOREACH()


SET(_output ${SCALAPACK_LIBRARY} ${BLACS_LIBRARIES})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCALAPACK DEFAULT_MSG
  _output # Cosmetic: Gives nice output
  SCALAPACK_LIBRARY
  LAPACK_FOUND
  )

IF(SCALAPACK_FOUND)
  SET(SCALAPACK_LIBRARIES
    ${SCALAPACK_LIBRARY}
    ${LAPACK_LIBRARIES}
    ${BLACS_LIBRARIES}
    )
  SET(SCALAPACK_LINKER_FLAGS
    ${LAPACK_LINKER_FLAGS}
    )

  MARK_AS_ADVANCED(
    lapack_LIBRARY
    atlas_LIBRARY
    blas_LIBRARY
    SCALAPACK_DIR
    SCALAPACK_LIBRARY
    BLACS_DIR
    BLACS_LIBRARY
    BLACSCINIT_LIBRARY
    BLACSF77INIT_LIBRARY
    )
ELSE()
  SET(SCALAPACK_DIR "" CACHE PATH
    "An optional hint to a SCALAPACK directory"
    )
  SET(BLACS_DIR "" CACHE PATH
    "An optional hint to a BLACS directory"
    )
ENDIF()

