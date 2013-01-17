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
# Try to find the petsc library
#
# This module exports:
#
#     PETSC_FOUND
#     PETSC_LIBRARIES
#     PETSC_INCLUDE_DIRS
#     PETSC_VERSION
#     PETSC_VERSION_MAJOR
#     PETSC_VERSION_MINOR
#     PETSC_VERSION_SUBMINOR
#     PETSC_VERSION_PATCH
#     PETSC_WITH_MPIUNI
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(PETSC_DIR "$ENV{PETSC_DIR}")
SET_IF_EMPTY(PETSC_ARCH "$ENV{PETSC_ARCH}")

#
# So, well, yes. I'd like to include the PETScConfig.cmake file via
# FIND_PACKAGE(), but it is broken beyond belief:
#
# - In source, i.e. PETSC_DIR/PETSC_ARCH, it sets BUILD_SHARED_LIBS.
# - It does not contain its very own version number
# - It does not contain its very own library location(s) or name(s)
# - It does not contain necessary includes
#
# - It writes a lot of FIND_LIBRARY(..) statements. Seriously. What the
#   heck? If its not the same library you're linking against, you cannot
#   assume to be API compatible, so why not just give a list of libraries?
#
# - It is not even considered to be installed in the gentoo petsc package
#   :-]
#

#
# TODO: We'll have to guess which external libraries we'll have to link
# against someday to avoid underlinkage
#

FIND_LIBRARY(PETSC_LIBRARY
  NAMES petsc
  HINTS
    # petsc is special. Account for that
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)


#
# So, up to this point it was easy. Now, the tricky part:
#

#
# Search for the first part of the includes:
#
FIND_PATH(PETSC_INCLUDE_DIR_ARCH petscconf.h
  HINTS
    # petsc is special. Account for that
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}
    ${PETSC_INCLUDE_DIRS}
  PATH_SUFFIXES petsc include include/petsc
)

#
# Sometimes, this is not enough...
# If petsc is not installed but in source tree layout, there will be
#   ${PETSC_DIR}/${PETSC_ARCH}/include - which we should have found by now.
#   ${PETSC_DIR}/include               - which we still have to find.
#
# Or it is installed in a non standard layout in the system (e.g. in
# Gentoo), where there will be
#   ${PETSC_DIR}/${PETSC_ARCH}/include
#   /usr/include/petsc ...
#
# Either way, we must be able to find petscversion.h:
#
FIND_PATH(PETSC_INCLUDE_DIR_COMMON petscversion.h
  HINTS
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}
    ${PETSC_INCLUDE_DIRS}
  PATH_SUFFIXES petsc include include/petsc
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(PETSC DEFAULT_MSG
  PETSC_LIBRARY
  PETSC_INCLUDE_DIR_ARCH
  PETSC_INCLUDE_DIR_COMMON
  )

IF(PETSC_FOUND)
  SET(PETSC_LIBRARIES
    ${PETSC_LIBRARY}
    )
  SET(PETSC_INCLUDE_DIRS
    ${PETSC_INCLUDE_DIR_ARCH}
    ${PETSC_INCLUDE_DIR_COMMON}
    )

  SET(PETSC_PETSCCONF_H "${PETSC_INCLUDE_DIR_ARCH}/petscconf.h")
  SET(PETSC_PETSCVERSION_H "${PETSC_INCLUDE_DIR_COMMON}/petscversion.h")

  #
  # Is petsc compiled with support for MPIUNI?
  #
  FILE(STRINGS "${PETSC_PETSCCONF_H}" PETSC_MPIUNI_STRING
    REGEX "#define.*PETSC_HAVE_MPIUNI 1")
  IF("${PETSC_MPIUNI_STRING}" STREQUAL "")
    SET(PETSC_WITH_MPIUNI FALSE)
  ELSE()
    SET(PETSC_WITH_MPIUNI TRUE)
  ENDIF()

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_MAJOR_STRING
    REGEX "#define.*PETSC_VERSION_MAJOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_MAJOR.*([0-9]+).*" "\\1"
    PETSC_VERSION_MAJOR "${PETSC_VERSION_MAJOR_STRING}"
    )

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_MINOR_STRING
    REGEX "#define.*PETSC_VERSION_MINOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_MINOR.*([0-9]+).*" "\\1"
    PETSC_VERSION_MINOR "${PETSC_VERSION_MINOR_STRING}"
    )

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_SUBMINOR_STRING
    REGEX "#define.*PETSC_VERSION_SUBMINOR")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_SUBMINOR.*([0-9]+).*" "\\1"
    PETSC_VERSION_SUBMINOR "${PETSC_VERSION_SUBMINOR_STRING}"
    )

  FILE(STRINGS "${PETSC_PETSCVERSION_H}" PETSC_VERSION_PATCH_STRING
    REGEX "#define.*PETSC_VERSION_PATCH")
  STRING(REGEX REPLACE "^.*PETSC_VERSION_PATCH.*([0-9]+).*" "\\1"
    PETSC_VERSION_PATCH "${PETSC_VERSION_PATCH_STRING}"
    )

  SET(PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}")

  MARK_AS_ADVANCED(
    PETSC_ARCH
    PETSC_DIR
    PETSC_INCLUDE_DIR_ARCH
    PETSC_INCLUDE_DIR_COMMON
    PETSC_INCLUDE_DIRS
    PETSC_LIBRARIES
  )
ELSE()
  SET(PETSC_DIR "" CACHE STRING
    "An optional hint to a PETSc directory"
    )
  SET(PETSC_ARCH "" CACHE STRING
    "An optional hint to a PETSc arch"
    )
ENDIF()

