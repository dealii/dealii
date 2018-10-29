## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
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
# Try to find the OpenCASCADE (OCC) library. This scripts supports the
# OpenCASCADE Community Edition (OCE) library, which is a cmake based
# OCC library. You might try the original OpenCASCADE library, but your
# mileage may vary.
#
# This module exports:
#
#   OPENCASCADE_DIR
#   OPENCASCADE_INCLUDE_DIRS
#   OPENCASCADE_LIBRARIES
#   OPENCASCADE_VERSION
#


SET(OPENCASCADE_DIR "" CACHE PATH "An optional hint to a OpenCASCADE installation")
SET_IF_EMPTY(OPENCASCADE_DIR "$ENV{OPENCASCADE_DIR}")
SET_IF_EMPTY(OPENCASCADE_DIR "$ENV{OCC_DIR}")
SET_IF_EMPTY(OPENCASCADE_DIR "$ENV{OCE_DIR}")
SET_IF_EMPTY(OPENCASCADE_DIR "$ENV{CASROOT}")


DEAL_II_FIND_PATH(OPENCASCADE_INCLUDE_DIR Standard_Version.hxx
  HINTS ${OPENCASCADE_DIR}
  PATH_SUFFIXES include include/oce include/opencascade inc
  )

IF(EXISTS ${OPENCASCADE_INCLUDE_DIR}/Standard_Version.hxx)
  FILE(STRINGS "${OPENCASCADE_INCLUDE_DIR}/Standard_Version.hxx" OPENCASCADE_VERSION
    REGEX "#define OCC_VERSION_COMPLETE "
    )
  STRING(REGEX REPLACE
    "#define OCC_VERSION_COMPLETE.*\"(.*)\"" "\\1"
    OPENCASCADE_VERSION "${OPENCASCADE_VERSION}"
    )
ENDIF()

# These seem to be pretty much the only required ones.
SET(_opencascade_libraries
  TKBO TKBool TKBRep TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo
  TKGeomBase TKHLR TKIGES TKMath TKMesh TKOffset TKPrim TKShHealing TKSTEP
  TKSTEPAttr TKSTEPBase TKSTEP209 TKSTL TKTopAlgo TKXSBase
  )

SET(_libraries "")
FOREACH(_library ${_opencascade_libraries})
  LIST(APPEND _libraries OPENCASCADE_${_library})
  DEAL_II_FIND_LIBRARY(OPENCASCADE_${_library}
    NAMES ${_library}
    HINTS ${OPENCASCADE_DIR}
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib mac64/clang/lib mac32/clang/lib lin64/gcc/lib lin32/gcc/lib
    )
ENDFOREACH()


DEAL_II_PACKAGE_HANDLE(OPENCASCADE
  LIBRARIES
    REQUIRED ${_libraries}
  INCLUDE_DIRS
    REQUIRED OPENCASCADE_INCLUDE_DIR
  USER_INCLUDE_DIRS
    REQUIRED OPENCASCADE_INCLUDE_DIR
  CLEAR
    _opencascade_libraries ${_libraries}
  )
