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
# Try to find the UMFPACK library
#
# This module exports
#
#   UMFPACK_LIBRARIES
#   UMFPACK_INCLUDE_DIRS
#   UMFPACK_LINKER_FLAGS
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(UMFPACK_DIR "$ENV{UMFPACK_DIR}")

#
# UMFPACK depends on AMD and BLAS, so search for them:
#
FIND_PACKAGE(AMD)
FIND_PACKAGE(BLAS)

FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h
  HINTS
    ${UMFPACK_DIR}
  PATH_SUFFIXES
    umfpack include/umfpack include Include UMFPACK/Include
  )

FIND_LIBRARY(UMFPACK_LIBRARY
  NAMES umfpack
  HINTS
    ${UMFPACK_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib Lib UMFPACK/Lib
  )

SET(required_variables
  BLAS_FOUND
  AMD_FOUND
  UMFPACK_INCLUDE_DIRS
  UMFPACK_LIBRARY
  )


#
# Well, recent versions of UMFPACK >= 5.6 include SuiteSparse_config.h, if so,
# ensure that we'll find these headers as well.
#
IF(NOT UMFPACK_INCLUDE_DIR MATCHES "-NOTFOUND")
  FILE(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_SUITESPARSE_STRING
    REGEX "#include \"SuiteSparse_config.h\"")
  IF(NOT "${UMFPACK_SUITESPARSE_STRING}" STREQUAL "")
    FIND_PACKAGE(SUITESPARSECONFIG)
    LIST(APPEND required_variables
      SUITESPARSECONFIG_FOUND
      )
  ENDIF()
ENDIF()
#
# Otherwise, we're lazy for the moment.
#


#
# UMFPACK may need to be linked with cholmod as well.
#
FIND_LIBRARY(CHOLMOD_LIBRARY
  NAMES cholmod
  HINTS
    ${UMFPACK_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib Lib UMFPACK/Lib
  )
IF(CHOLMOD_LIBRARY MATCHES "-NOTFOUND")
  SET(CHOLMOD_LIBRARY "")
ENDIF()


SET(UMFPACK_LIBRARIES
  ${UMFPACK_LIBRARY}
  ${CHOLMOD_LIBRARY} # may be empty
  ${AMD_LIBRARY}
  ${SUITESPARSECONFIG_LIBRARY} # may be empty
  ${BLAS_LIBRARIES}
  )

SET(UMFPACK_INCLUDE_DIRS
  ${UMFPACK_INCLUDE_DIR}
  ${SUITESPARSECONFIG_INCLUDE_DIR} # may be empty
  ${AMD_INCLUDE_DIR}
  )

SET(UMFPACK_LINKER_FLAGS
  ${BLAS_LINKER_FLAGS}
  )


FIND_PACKAGE_HANDLE_STANDARD_ARGS(UMFPACK DEFAULT_MSG ${required_variables})

IF(UMFPACK_FOUND)
  MARK_AS_ADVANCED(
    UMFPACK_LIBRARY
    UMFPACK_INCLUDE_DIR
    UMFPACK_DIR
    CHOLMOD_LIBRARY
    atlas_LIBRARY
    blas_LIBRARY
  )
ELSE()
  SET(UMFPACK_DIR "" CACHE STRING
    "An optional hint to an UMFPACK directory"
    )
ENDIF()

