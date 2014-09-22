## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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

SET(SLEPC_DIR "" CACHE PATH "An optional hint to a SLEPC directory")
SET_IF_EMPTY(SLEPC_DIR "$ENV{SLEPC_DIR}")
SET_IF_EMPTY(PETSC_DIR "$ENV{PETSC_DIR}")
SET_IF_EMPTY(PETSC_ARCH "$ENV{PETSC_ARCH}")

#
# Luckily, SLEPc wants the same insanity as PETSc, so we can just copy the
# mechanism.
#

DEAL_II_FIND_LIBRARY(SLEPC_LIBRARY
  NAMES slepc
  HINTS ${SLEPC_DIR} ${SLEPC_DIR}/${PETSC_ARCH} ${PETSC_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)

DEAL_II_FIND_PATH(SLEPC_INCLUDE_DIR_ARCH slepcconf.h
  HINTS
    ${SLEPC_DIR}
    ${SLEPC_DIR}/${PETSC_ARCH}
    ${SLEPC_INCLUDE_DIRS}
    ${PETSC_DIR}
  PATH_SUFFIXES slepc include include/slepc
)

DEAL_II_FIND_PATH(SLEPC_INCLUDE_DIR_COMMON slepcversion.h
  HINTS
    ${SLEPC_DIR}
    ${SLEPC_DIR}/${PETSC_ARCH}
    ${SLEPC_INCLUDE_DIRS}
    ${PETSC_DIR}
  PATH_SUFFIXES slepc include include/slepc
)

SET(SLEPC_SLEPCVERSION_H "${SLEPC_INCLUDE_DIR_COMMON}/slepcversion.h")
IF(EXISTS ${SLEPC_SLEPCVERSION_H})
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
  SET(SLEPC_VERSION
    "${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}.${SLEPC_VERSION_SUBMINOR}.${SLEPC_VERSION_PATCH}"
    )
ENDIF()

DEAL_II_PACKAGE_HANDLE(SLEPC
  LIBRARIES
    REQUIRED SLEPC_LIBRARY PETSC_LIBRARIES
  INCLUDE_DIRS
    REQUIRED SLEPC_INCLUDE_DIR_ARCH SLEPC_INCLUDE_DIR_COMMON
  USER_INCLUDE_DIRS
    REQUIRED SLEPC_INCLUDE_DIR_ARCH SLEPC_INCLUDE_DIR_COMMON
  CLEAR SLEPC_LIBRARY SLEPC_INCLUDE_DIR_ARCH SLEPC_INCLUDE_DIR_COMMON
  )
