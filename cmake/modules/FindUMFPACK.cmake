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
# Try to find the UMFPACK library
#
# This module exports
#
#   UMFPACK_LIBRARIES
#   UMFPACK_INCLUDE_DIRS
#   UMFPACK_LINKER_FLAGS
#   UMFPACK_VERSION
#   UMFPACK_VERSION_MAJOR
#   UMFPACK_VERSION_MINOR
#   UMFPACK_VERSION_SUBMINOR
#

SET(UMFPACK_DIR "" CACHE PATH "An optional hint to an UMFPACK directory")
SET(SUITESPARSE_DIR "" CACHE PATH
  "An optional hint to a SUITESPARSE directory"
  )
FOREACH(_comp SUITESPARSE SUITESPARSE_CONFIG UMFPACK AMD CHOLMOD COLAMD)
  SET_IF_EMPTY(${_comp}_DIR "$ENV{${_comp}_DIR}")
ENDFOREACH()

#
# Two macros to make life easier:
#

MACRO(FIND_UMFPACK_PATH _comp _file)
  STRING(TOLOWER ${_comp} _comp_lowercase)
  STRING(TOUPPER ${_comp} _comp_uppercase)
  DEAL_II_FIND_PATH(${_comp}_INCLUDE_DIR ${_file}
    HINTS
      ${${_comp_uppercase}_DIR}
      ${SUITESPARSE_DIR}/${_comp}
      ${UMFPACK_DIR}/../${_comp}
      ${UMFPACK_DIR}/${_comp}
      ${UMFPACK_DIR}
    PATH_SUFFIXES
      ${_comp_lowercase} include/${_comp_lowercase} include Include ${_comp}/Include suitesparse
    )
ENDMACRO()

MACRO(FIND_UMFPACK_LIBRARY _comp _name)
  STRING(TOUPPER ${_comp} _comp_uppercase)
  DEAL_II_FIND_LIBRARY(${_comp}_LIBRARY
    NAMES ${_name}
    HINTS
      ${${_comp_uppercase}_DIR}
      ${SUITESPARSE_DIR}/${_comp}
      ${UMFPACK_DIR}/../${_comp}
      ${UMFPACK_DIR}/${_comp}
      ${UMFPACK_DIR}
    PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib Lib ${_comp}/Lib
    )
ENDMACRO()


#
# Search for include directories:
#
FIND_UMFPACK_PATH(UMFPACK umfpack.h)
FIND_UMFPACK_PATH(AMD amd.h)

IF(EXISTS ${UMFPACK_INCLUDE_DIR}/umfpack.h)
  #
  # Well, recent versions of UMFPACK include SuiteSparse_config.h, if so,
  # ensure that we'll find these headers as well.
  #
  FILE(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_SUITESPARSE_STRING
    REGEX "#include \"SuiteSparse_config.h\"")
  IF(NOT "${UMFPACK_SUITESPARSE_STRING}" STREQUAL "")
    FIND_UMFPACK_PATH(SuiteSparse_config SuiteSparse_config.h)
  ENDIF()

  FILE(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_VERSION_MAJOR_STRING
    REGEX "#define.*UMFPACK_MAIN_VERSION")
  STRING(REGEX REPLACE "^.*UMFPACK_MAIN_VERSION.*([0-9]+).*" "\\1"
    UMFPACK_VERSION_MAJOR "${UMFPACK_VERSION_MAJOR_STRING}"
    )
  FILE(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_VERSION_MINOR_STRING
    REGEX "#define.*UMFPACK_SUB_VERSION")
  STRING(REGEX REPLACE "^.*UMFPACK_SUB_VERSION.*([0-9]+).*" "\\1"
    UMFPACK_VERSION_MINOR "${UMFPACK_VERSION_MINOR_STRING}"
    )
  FILE(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_VERSION_SUBMINOR_STRING
    REGEX "#define.*UMFPACK_SUBSUB_VERSION")
  STRING(REGEX REPLACE "^.*UMFPACK_SUBSUB_VERSION.*([0-9]+).*" "\\1"
    UMFPACK_VERSION_SUBMINOR "${UMFPACK_VERSION_SUBMINOR_STRING}"
    )
  SET(UMFPACK_VERSION
    "${UMFPACK_VERSION_MAJOR}.${UMFPACK_VERSION_MINOR}.${UMFPACK_VERSION_SUBMINOR}"
    )
ENDIF()

#
# Link against everything we can find to avoid underlinkage:
#
FIND_UMFPACK_LIBRARY(UMFPACK umfpack)
FIND_UMFPACK_LIBRARY(AMD amd)
FIND_UMFPACK_LIBRARY(CHOLMOD cholmod)
FIND_UMFPACK_LIBRARY(COLAMD colamd)
FIND_UMFPACK_LIBRARY(CCOLAMD ccolamd)
FIND_UMFPACK_LIBRARY(CAMD camd)
FIND_UMFPACK_LIBRARY(SuiteSparse_config suitesparseconfig)

#
# Test whether libsuitesparseconfig.xxx can be used for shared library
# linkage. If not, exclude it from the command line.
#
LIST(APPEND CMAKE_REQUIRED_LIBRARIES
  "-shared"
  ${SuiteSparse_config_LIBRARY}
  )
CHECK_CXX_SOURCE_COMPILES("extern int SuiteSparse_version (int[3]);
  void foo(int bar[3]) { SuiteSparse_version(bar);}"
  LAPACK_SUITESPARSECONFIG_WITH_PIC
  )
RESET_CMAKE_REQUIRED()

IF(LAPACK_SUITESPARSECONFIG_WITH_PIC OR NOT BUILD_SHARED_LIBS)
  SET(_suitesparse_config SuiteSparse_config_LIBRARY)
ENDIF()

#
# Add rt to the link interface as well (for whatever reason,
# libsuitesparse.so depends on clock_gettime but the shared
# lib does not record its dependence on librt.so as evidenced
# by ldd :-( ):
#
FIND_SYSTEM_LIBRARY(rt_LIBRARY NAMES rt)
MARK_AS_ADVANCED(rt_LIBRARY)

DEAL_II_PACKAGE_HANDLE(UMFPACK
  LIBRARIES
    REQUIRED UMFPACK_LIBRARY
    OPTIONAL CHOLMOD_LIBRARY CCOLAMD_LIBRARY COLAMD_LIBRARY CAMD_LIBRARY ${_suitesparse_config}
    REQUIRED AMD_LIBRARY
    OPTIONAL METIS_LIBRARIES LAPACK_LIBRARIES rt_LIBRARY
  INCLUDE_DIRS
    REQUIRED UMFPACK_INCLUDE_DIR AMD_INCLUDE_DIR
    OPTIONAL SuiteSparse_config_INCLUDE_DIR
  LINKER_FLAGS
    OPTIONAL LAPACK_LINKER_FLAGS
  CLEAR
    UMFPACK_LIBRARY CHOLMOD_LIBRARY CCOLAMD_LIBRARY COLAMD_LIBRARY
    CAMD_LIBRARY SuiteSparse_config_LIBRARY AMD_LIBRARY UMFPACK_INCLUDE_DIR
    AMD_INCLUDE_DIR SuiteSparse_config_INCLUDE_DIR
  )

IF(UMFPACK_FOUND)
  MARK_AS_ADVANCED(SUITESPARSE_DIR)
ELSE()
  MARK_AS_ADVANCED(CLEAR SUITESPARSE_DIR)
ENDIF()
