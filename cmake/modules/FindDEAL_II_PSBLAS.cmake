## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2025 by the deal.II authors
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
# Try to find the PSBLAS library
#
# This module exports
#
#   PSBLAS_LIBRARY
#   PSBLAS_INCLUDE_DIR
#

set(PSBLAS_DIR "" CACHE PATH "An optional hint to a PSBLAS installation containing the PSBLAS include directory and libraries")
set_if_empty(PSBLAS_DIR "$ENV{PSBLAS_DIR}")

set(_psblas_libs "psb_base;psb_cbind;psb_krylov;psb_prec;psb_util")
set(_psblas_library_variables "")

foreach(_lib ${_psblas_libs})
  string(TOUPPER ${_lib} _lib_upper)
  string(REPLACE "PSB_" "PSBLAS_" _var_name "${_lib_upper}_LIBRARY")
  deal_ii_find_library(${_var_name}
    NAMES ${_lib}
    HINTS ${PSBLAS_DIR}
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )
  list(APPEND _psblas_library_variables ${_var_name})
endforeach()

deal_ii_find_path(PSBLAS_INCLUDE_DIR psb_c_base.h
  HINTS ${PSBLAS_DIR}
  PATH_SUFFIXES include
  )

set(_additional_libraries "")
set(_fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
# Arbitrarily set this to gfortran and m, if they are not found, the user will
# have to set them manually
set_if_empty(_fortran_libs gfortran quadmath m)

foreach(_lib ${_fortran_libs})
  find_system_library(${_lib}_LIBRARY NAMES ${_lib})
  list(APPEND _additional_libraries ${_lib}_LIBRARY)
endforeach()

process_feature(PSBLAS
  LIBRARIES 
    REQUIRED 
      ${_psblas_library_variables}
      ${_additional_libraries}
  INCLUDE_DIRS 
    REQUIRED PSBLAS_INCLUDE_DIR
  CLEAR
    ${_psblas_library_variables} ${_additional_libraries} PSBLAS_INCLUDE_DIR
  )
