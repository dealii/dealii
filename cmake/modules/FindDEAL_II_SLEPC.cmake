## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
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

set(SLEPC_DIR "" CACHE PATH "An optional hint to a SLEPC directory")
set_if_empty(SLEPC_DIR "$ENV{SLEPC_DIR}")
set_if_empty(PETSC_DIR "$ENV{PETSC_DIR}")
set_if_empty(PETSC_ARCH "$ENV{PETSC_ARCH}")

#
# Luckily, SLEPc wants the same insanity as PETSc, so we can just copy the
# mechanism.
#

deal_ii_find_library(SLEPC_LIBRARY
  NAMES slepc libslepc
  HINTS ${SLEPC_DIR} ${SLEPC_DIR}/${PETSC_ARCH} ${PETSC_DIR} ${PETSC_DIR}/${PETSC_ARCH}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)

deal_ii_find_path(SLEPC_INCLUDE_DIR_ARCH slepcconf.h
  HINTS
    ${SLEPC_DIR}
    ${SLEPC_DIR}/${PETSC_ARCH}
    ${SLEPC_INCLUDE_DIRS}
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}
  PATH_SUFFIXES slepc include include/slepc
)

deal_ii_find_path(SLEPC_INCLUDE_DIR_COMMON slepcversion.h
  HINTS
    ${SLEPC_DIR}
    ${SLEPC_DIR}/${PETSC_ARCH}
    ${SLEPC_INCLUDE_DIRS}
    ${PETSC_DIR}
    ${PETSC_DIR}/${PETSC_ARCH}
  PATH_SUFFIXES slepc include include/slepc
)

set(SLEPC_SLEPCVERSION_H "${SLEPC_INCLUDE_DIR_COMMON}/slepcversion.h")
if(EXISTS ${SLEPC_SLEPCVERSION_H})
  file(STRINGS "${SLEPC_SLEPCVERSION_H}" SLEPC_VERSION_MAJOR_STRING
    REGEX "#define.*SLEPC_VERSION_MAJOR")
  string(REGEX REPLACE "^.*SLEPC_VERSION_MAJOR.* ([0-9]+).*" "\\1"
    SLEPC_VERSION_MAJOR "${SLEPC_VERSION_MAJOR_STRING}"
    )
  file(STRINGS "${SLEPC_SLEPCVERSION_H}" SLEPC_VERSION_MINOR_STRING
    REGEX "#define.*SLEPC_VERSION_MINOR")
  string(REGEX REPLACE "^.*SLEPC_VERSION_MINOR.* ([0-9]+).*" "\\1"
    SLEPC_VERSION_MINOR "${SLEPC_VERSION_MINOR_STRING}"
    )
  file(STRINGS "${SLEPC_SLEPCVERSION_H}" SLEPC_VERSION_SUBMINOR_STRING
    REGEX "#define.*SLEPC_VERSION_SUBMINOR")
  string(REGEX REPLACE "^.*SLEPC_VERSION_SUBMINOR.* ([0-9]+).*" "\\1"
    SLEPC_VERSION_SUBMINOR "${SLEPC_VERSION_SUBMINOR_STRING}"
    )
  file(STRINGS "${SLEPC_SLEPCVERSION_H}" SLEPC_VERSION_PATCH_STRING
    REGEX "#define.*SLEPC_VERSION_PATCH")
  string(REGEX REPLACE "^.*SLEPC_VERSION_PATCH.* ([0-9]+).*" "\\1"
    SLEPC_VERSION_PATCH "${SLEPC_VERSION_PATCH_STRING}"
    )
  set(SLEPC_VERSION
    "${SLEPC_VERSION_MAJOR}.${SLEPC_VERSION_MINOR}.${SLEPC_VERSION_SUBMINOR}.${SLEPC_VERSION_PATCH}"
    )
endif()

process_feature(SLEPC
  LIBRARIES
    REQUIRED SLEPC_LIBRARY PETSC_LIBRARIES
  INCLUDE_DIRS
    REQUIRED SLEPC_INCLUDE_DIR_ARCH SLEPC_INCLUDE_DIR_COMMON
  CLEAR SLEPC_LIBRARY SLEPC_INCLUDE_DIR_ARCH SLEPC_INCLUDE_DIR_COMMON
  )
