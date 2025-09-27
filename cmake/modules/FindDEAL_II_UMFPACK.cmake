## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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

set(UMFPACK_DIR "" CACHE PATH "An optional hint to an UMFPACK directory")
set(SUITESPARSE_DIR "" CACHE PATH
  "An optional hint to a SUITESPARSE directory"
  )
foreach(_comp SUITESPARSE SUITESPARSE_CONFIG UMFPACK AMD CHOLMOD COLAMD)
  set_if_empty(${_comp}_DIR "$ENV{${_comp}_DIR}")
endforeach()

#
# Two macros to make life easier:
#

macro(find_umfpack_path _comp _file)
  string(TOLOWER ${_comp} _comp_lowercase)
  string(TOUPPER ${_comp} _comp_uppercase)
  deal_ii_find_path(${_comp}_INCLUDE_DIR ${_file}
    HINTS
      ${${_comp_uppercase}_DIR}
      ${SUITESPARSE_DIR}/${_comp}
      ${UMFPACK_DIR}/../${_comp}
      ${UMFPACK_DIR}/${_comp}
      ${UMFPACK_DIR}
    PATH_SUFFIXES
      ${_comp_lowercase} include/${_comp_lowercase} include Include ${_comp}/Include suitesparse
    )
endmacro()

macro(find_umfpack_library _comp _name)
  string(TOUPPER ${_comp} _comp_uppercase)
  deal_ii_find_library(${_comp}_LIBRARY
    NAMES ${_name} lib${_name}
    HINTS
      ${${_comp_uppercase}_DIR}
      ${SUITESPARSE_DIR}/${_comp}
      ${UMFPACK_DIR}/../${_comp}
      ${UMFPACK_DIR}/${_comp}
      ${UMFPACK_DIR}
    PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib Lib ${_comp}/Lib
    )
endmacro()


#
# Search for include directories:
#
FIND_UMFPACK_PATH(UMFPACK umfpack.h)
FIND_UMFPACK_PATH(AMD amd.h)

if(EXISTS ${UMFPACK_INCLUDE_DIR}/umfpack.h)
  #
  # Well, recent versions of UMFPACK include SuiteSparse_config.h, if so,
  # ensure that we'll find these headers as well.
  #
  file(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_SUITESPARSE_STRING
    REGEX "#include \"SuiteSparse_config.h\"")
  if(NOT "${UMFPACK_SUITESPARSE_STRING}" STREQUAL "")
    FIND_UMFPACK_PATH(SuiteSparse_config SuiteSparse_config.h)
  endif()

  file(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_VERSION_MAJOR_STRING
    REGEX "#define.*UMFPACK_MAIN_VERSION")
  string(REGEX REPLACE "^.*UMFPACK_MAIN_VERSION.*([0-9]+).*" "\\1"
    UMFPACK_VERSION_MAJOR "${UMFPACK_VERSION_MAJOR_STRING}"
    )
  file(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_VERSION_MINOR_STRING
    REGEX "#define.*UMFPACK_SUB_VERSION")
  string(REGEX REPLACE "^.*UMFPACK_SUB_VERSION.*([0-9]+).*" "\\1"
    UMFPACK_VERSION_MINOR "${UMFPACK_VERSION_MINOR_STRING}"
    )
  file(STRINGS "${UMFPACK_INCLUDE_DIR}/umfpack.h" UMFPACK_VERSION_SUBMINOR_STRING
    REGEX "#define.*UMFPACK_SUBSUB_VERSION")
  string(REGEX REPLACE "^.*UMFPACK_SUBSUB_VERSION.*([0-9]+).*" "\\1"
    UMFPACK_VERSION_SUBMINOR "${UMFPACK_VERSION_SUBMINOR_STRING}"
    )
  set(UMFPACK_VERSION
    "${UMFPACK_VERSION_MAJOR}.${UMFPACK_VERSION_MINOR}.${UMFPACK_VERSION_SUBMINOR}"
    )
endif()

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
# Add rt to the link interface as well (for whatever reason,
# libsuitesparse.so depends on clock_gettime but the shared
# lib does not record its dependence on librt.so as evidenced
# by ldd :-( ):
#
find_system_library(rt_LIBRARY NAMES rt)

process_feature(UMFPACK
  LIBRARIES
    REQUIRED UMFPACK_LIBRARY
    OPTIONAL CHOLMOD_LIBRARY CCOLAMD_LIBRARY COLAMD_LIBRARY CAMD_LIBRARY SuiteSparse_config_LIBRARY
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

if(UMFPACK_FOUND)
  mark_as_advanced(SUITESPARSE_DIR)
else()
  mark_as_advanced(CLEAR SUITESPARSE_DIR)
endif()
