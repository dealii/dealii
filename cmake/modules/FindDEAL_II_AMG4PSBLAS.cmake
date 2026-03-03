## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
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
# Try to find the AMG4PSBLAS library
#
# This module exports
#
#   
#   AMG4PSBLAS_INCLUDE_DIR
#

set(AMG4PSBLAS_DIR "" CACHE PATH "An optional hint to a AMG4PSBLAS installation containing the AMG4PSBLAS include directory and libraries")
set_if_empty(AMG4PSBLAS_DIR "$ENV{AMG4PSBLAS_DIR}")

set(_amg4psblas_libs "amg_cbind;amg_prec")
set(_amg4psblas_library_variables "")

foreach(_lib ${_amg4psblas_libs})
  string(TOUPPER ${_lib} _lib_upper)
  string(REPLACE "AMG_" "AMG4PSBLAS_" _var_name "${_lib_upper}_LIBRARY")
  deal_ii_find_library(${_var_name}
    NAMES ${_lib}
    HINTS ${AMG4PSBLAS_DIR}
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )
  list(APPEND _amg4psblas_library_variables ${_var_name})
endforeach()

deal_ii_find_path(AMG4PSBLAS_INCLUDE_DIR amg_cbind.h
  HINTS ${AMG4PSBLAS_DIR}
  PATH_SUFFIXES include
  )

  message(STATUS "AMG4PSBLAS_INCLUDE_DIR: ${AMG4PSBLAS_INCLUDE_DIR}")


set(AMG4PSBLAS_PSBLASVERSION_H "${AMG4PSBLAS_INCLUDE_DIR}/amg_config.h")
if(EXISTS ${AMG4PSBLAS_PSBLASVERSION_H})
  file(STRINGS "${AMG4PSBLAS_PSBLASVERSION_H}" AMG4PSBLAS_VERSION_MAJOR_STRING
    REGEX "^#[ \t]*define[ \t]+AMG_VERSION_MAJOR[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+AMG_VERSION_MAJOR[ \t]+([0-9]+)[ \t]*$" "\\1"
    AMG4PSBLAS_VERSION_MAJOR "${AMG4PSBLAS_VERSION_MAJOR_STRING}")

  file(STRINGS "${AMG4PSBLAS_PSBLASVERSION_H}" AMG4PSBLAS_VERSION_MINOR_STRING
    REGEX "^#[ \t]*define[ \t]+AMG_VERSION_MINOR[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+AMG_VERSION_MINOR[ \t]+([0-9]+)[ \t]*$" "\\1"
    AMG4PSBLAS_VERSION_MINOR "${AMG4PSBLAS_VERSION_MINOR_STRING}")

  file(STRINGS "${AMG4PSBLAS_PSBLASVERSION_H}" AMG4PSBLAS_VERSION_PATCHLEVEL_STRING
    REGEX "^#[ \t]*define[ \t]+AMG_VERSION_PATCHLEVEL[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+AMG_VERSION_PATCHLEVEL[ \t]+([0-9]+)[ \t]*$" "\\1"
    AMG4PSBLAS_VERSION_PATCHLEVEL "${AMG4PSBLAS_VERSION_PATCHLEVEL_STRING}")

  set(AMG4PSBLAS_VERSION
    "${AMG4PSBLAS_VERSION_MAJOR}.${AMG4PSBLAS_VERSION_MINOR}.${AMG4PSBLAS_VERSION_PATCHLEVEL}"
    )
  message(STATUS "AMG4PSBLAS version detected: ${AMG4PSBLAS_VERSION}")
endif()

set(_additional_libraries ${_interface_lapack} ${_interface_blas} ${_interface_psblas})
set(_fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
# Arbitrarily set this to gfortran and m, if they are not found, the user will
# have to set them manually
set_if_empty(_fortran_libs gfortran quadmath m)

foreach(_lib ${_fortran_libs})
  find_system_library(${_lib}_LIBRARY NAMES ${_lib})
  list(APPEND _additional_libraries ${_lib}_LIBRARY)
endforeach()

process_feature(AMG4PSBLAS
  LIBRARIES 
    REQUIRED 
    ${_amg4psblas_library_variables}
    ${_additional_libraries}
    LAPACK_LIBRARIES
    OPTIONAL
    MPI_CXX_LIBRARIES
    MPI_Fortran_LIBRARIES
  INCLUDE_DIRS 
    REQUIRED AMG4PSBLAS_INCLUDE_DIR
  LINKER_FLAGS
    REQUIRED ${_interface_lapack} ${_interface_blas}
  CLEAR
    ${_amg4psblas_library_variables} ${_additional_libraries} AMG4PSBLAS_INCLUDE_DIR
  )