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
# Try to find the SLEPC library
#
# This module exports:
#
#     SLEPC_FOUND
#     SLEPC_LIBRARIES
#     SLEPC_INCLUDE_DIRS
#     SLEPC_VERSION
#     SLEPC_VERSION_MAJOR
#     SLEPC_VERSION_MINOR
#     SLEPC_VERSION_SUBMINOR
#     SLEPC_VERSION_PATCH
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(SLEPC_DIR "$ENV{SLEPC_DIR}")
SET_IF_EMPTY(SLEPC_ARCH "$ENV{SLEPC_ARCH}")
SET_IF_EMPTY(PETSC_DIR "$ENV{PETSC_DIR}")

#
# Luckily, SLEPc wants the same insanity as PETSc, so we can just copy the
# mechanism.
#

FIND_LIBRARY(SLEPC_LIBRARIES
  NAMES slepc
  HINTS
    # SLEPC is special. Account for that
    ${SLEPC_DIR}
    ${SLEPC_DIR}/${SLEPC_ARCH}
    ${PETSC_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)


#
# So, up to this point it was easy. Now, the tricky part:
#


#
# Search for the first part of the includes:
#
FIND_PATH(SLEPC_INCLUDE_DIR_ARCH slepcconf.h
  HINTS
    # SLEPC is special. Account for that
    ${SLEPC_DIR}
    ${SLEPC_DIR}/${SLEPC_ARCH}
    ${SLEPC_INCLUDE_DIRS}
    ${PETSC_DIR}
  PATH_SUFFIXES slepc include include/slepc
)

#
# Sometimes, this is not enough...
# If SLEPc is not installed but in source tree layout, there will be
#   ${SLEPC_DIR}/${SLEPC_ARCH}/include - which we should have found by now.
#   ${SLEPC_DIR}/include               - which we still have to find.
#
# Or it is installed in a non standard layout in the system (e.g. in
# Gentoo), where there will be
#   ${SLEPC_DIR}/${SLEPC_ARCH}/include
#   /usr/include/slepc ...
#
# Either way, slepcversion.h should lie around:
#
FIND_PATH(SLEPC_INCLUDE_DIR_COMMON slepcversion.h
  HINTS
    ${SLEPC_DIR}
    ${SLEPC_DIR}/${SLEPC_ARCH}
    ${SLEPC_INCLUDE_DIRS}
    ${PETSC_DIR}
  PATH_SUFFIXES slepc include include/slepc
)

#
# And finally set SLEPC_INCLUDE_DIRS depending on the outcome of our crude
# guess:
#
IF( SLEPC_INCLUDE_DIR_ARCH MATCHES "-NOTFOUND" OR
    SLEPC_INCLUDE_DIR_COMMON MATCHES "-NOTFOUND" )
  SET(SLEPC_INCLUDE_DIRS "SLEPC_INCLUDE_DIRS-NOTFOUND"
    CACHE STRING "Include paths for SLEPC"
    FORCE
    )
  UNSET(SLEPC_INCLUDE_DIR_ARCH CACHE)
  UNSET(SLEPC_INCLUDE_DIR_COMMON CACHE)
ELSE()
  UNSET(SLEPC_INCLUDE_DIRS CACHE)
  SET(SLEPC_INCLUDE_DIRS
    ${SLEPC_INCLUDE_DIR_ARCH}
    ${SLEPC_INCLUDE_DIR_COMMON}
    )
ENDIF()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SLEPC DEFAULT_MSG
  SLEPC_LIBRARIES
  SLEPC_INCLUDE_DIRS
  )

IF(SLEPC_FOUND)
  SET(SLEPC_SLEPCVERSION_H "${SLEPC_INCLUDE_DIR_COMMON}/slepcversion.h")

  FILE(STRINGS "${SLEPC_SLEPCVERSION_H}" SLEPC_VERSION_MAJOR_STRING
    REGEX "#define.*SLEPC_VERSION_MAJOR")
  STRING(REGEX REPLACE "^.*SLEPC_VERSION_MAJOR.*([0-9]+).*" "\\1"
    SLEPC_VERSION_MAJOR "${SLEPC_VERSION_MAJOR_STRING}"
    )

  FILE(STRINGS "${SLEPC_SLEPCVERSION_H}" SLEPC_VERSION_MINOR_STRING
    REGEX "#define.*SLEPC_VERSION_MINOR")
  STRING(REGEX REPLACE "^.*SLEPC_VERSION_MINOR.*([0-9]+).*" "\\1"
    SLEPC_VERSION_MINOR "${SLEPC_VERSION_MINOR_STRING}"
    )

  FILE(STRINGS "${SLEPC_SLEPCVERSION_H}" SLEPC_VERSION_SUBMINOR_STRING
    REGEX "#define.*SLEPC_VERSION_SUBMINOR")
  STRING(REGEX REPLACE "^.*SLEPC_VERSION_SUBMINOR.*([0-9]+).*" "\\1"
    SLEPC_VERSION_SUBMINOR "${SLEPC_VERSION_SUBMINOR_STRING}"
    )

  FILE(STRINGS "${SLEPC_SLEPCVERSION_H}" SLEPC_VERSION_PATCH_STRING
    REGEX "#define.*SLEPC_VERSION_PATCH")
  STRING(REGEX REPLACE "^.*SLEPC_VERSION_PATCH.*([0-9]+).*" "\\1"
    SLEPC_VERSION_PATCH "${SLEPC_VERSION_PATCH_STRING}"
    )

  SET(SLEPC_VERSION "${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}.${SLEPC_VERSION_SUBMINOR}")

  MARK_AS_ADVANCED(
    SLEPC_ARCH
    SLEPC_DIR
    SLEPC_INCLUDE_DIR_ARCH
    SLEPC_INCLUDE_DIR_COMMON
    SLEPC_INCLUDE_DIRS
    SLEPC_LIBRARIES
  )
ELSE()
  SET(SLEPC_DIR "" CACHE STRING
    "An optional hint to a SLEPC directory"
    )
  SET(SLEPC_ARCH "" CACHE STRING
    "An optional hint to a SLEPC arch"
    )
ENDIF()

