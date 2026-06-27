## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2017 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Try to find the PSBLAS library
#
# This module exports
#
#   PSBLAS_LIBRARY
#   PSBLAS_INCLUDE_DIR
#   PSBLAS_VERSION
#   PSB_VERSION_MAJOR
#   PSB_VERSION_MINOR
#   PSB_VERSION_PATCHLEVEL
#

set(PSBLAS_DIR "" CACHE PATH "An optional hint to a PSBLAS installation containing the PSBLAS include directory and libraries")
set_if_empty(PSBLAS_DIR "$ENV{PSBLAS_DIR}")

set(_psblas_libs "psb_base;psb_cbind;psb_linsolve;psb_prec;psb_ext;psb_util")
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

set(PSBLAS_PSBLASVERSION_H "${PSBLAS_INCLUDE_DIR}/psb_config.h")
if(EXISTS ${PSBLAS_PSBLASVERSION_H})
  file(STRINGS "${PSBLAS_PSBLASVERSION_H}" PSB_VERSION_MAJOR_STRING
    REGEX "^#[ \t]*define[ \t]+PSB_VERSION_MAJOR[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+PSB_VERSION_MAJOR[ \t]+([0-9]+)[ \t]*$" "\\1"
    PSB_VERSION_MAJOR "${PSB_VERSION_MAJOR_STRING}"
    )
  file(STRINGS "${PSBLAS_PSBLASVERSION_H}" PSB_VERSION_MINOR_STRING
    REGEX "^#[ \t]*define[ \t]+PSB_VERSION_MINOR[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+PSB_VERSION_MINOR[ \t]+([0-9]+)[ \t]*$" "\\1"
    PSB_VERSION_MINOR "${PSB_VERSION_MINOR_STRING}"
    )
  file(STRINGS "${PSBLAS_PSBLASVERSION_H}" PSB_VERSION_PATCHLEVEL_STRING
    REGEX "^#[ \t]*define[ \t]+PSB_VERSION_PATCHLEVEL[ \t]+[0-9]+[ \t]*$")
  string(REGEX REPLACE "^#[ \t]*define[ \t]+PSB_VERSION_PATCHLEVEL[ \t]+([0-9]+)[ \t]*$" "\\1"
    PSB_VERSION_PATCHLEVEL "${PSB_VERSION_PATCHLEVEL_STRING}"
    )
  set(PSBLAS_VERSION
    "${PSB_VERSION_MAJOR}.${PSB_VERSION_MINOR}.${PSB_VERSION_PATCHLEVEL}"
    )
endif()

#
# PSBLAS is a Fortran library distributed as static archives. We have to
# manually pick up the complete link interface. If
# CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is not available, do it
# unconditionally for the most common case (gfortran).
#
# Since CMake 3.9 the gcc runtime library libgcc.a has been added to the
# CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES variable. Additionally macOS GCC
# adds emutls_w and heapt_w. We remove all of these because explicitly
# linking them from a C++ linker invocation can break the build.
#
set(_fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
list(REMOVE_ITEM _fortran_libs gcc_s gcc emutls_w heapt_w)
set_if_empty(_fortran_libs gfortran quadmath m)

set(_additional_libraries "")
# Platform-specific shared/static library extensions for the
# compiler -print-file-name lookup.
if(APPLE)
  set(_psblas_lib_extensions ".so" ".so.1" ".dylib" ".a")
elseif(WIN32)
  set(_psblas_lib_extensions ".dll" ".lib" ".a")
else()
  set(_psblas_lib_extensions ".so" ".so.1" ".a" ".dylib")
endif()

foreach(_lib ${_fortran_libs})
  #
  # find_system_library uses check_cxx_source_compiles with "-l<name>" which
  # fails when the C++ compiler wrapper does not know the Fortran runtime's
  # library directory. We therefore also try to resolve the library via the
  # Fortran compiler's -print-file-name.
  #
  find_system_library(${_lib}_LIBRARY NAMES ${_lib})
  if(NOT ${_lib}_LIBRARY OR ${_lib}_LIBRARY MATCHES "NOTFOUND")
    foreach(__ext ${_psblas_lib_extensions})
      execute_process(
        COMMAND ${CMAKE_Fortran_COMPILER}
                -print-file-name=lib${_lib}${__ext}
        OUTPUT_VARIABLE __candidate
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        )
      if(__candidate AND EXISTS "${__candidate}")
        set(${_lib}_LIBRARY "${__candidate}" CACHE FILEPATH
            "Location of ${_lib} runtime library" FORCE)
        set(${_lib}_LIBRARY "${__candidate}")
        mark_as_advanced(${_lib}_LIBRARY)
        break()
      endif()
    endforeach()
  endif()
  list(APPEND _additional_libraries ${_lib}_LIBRARY)
endforeach()

process_feature(PSBLAS
  LIBRARIES
    REQUIRED
      ${_psblas_library_variables}
      ${_additional_libraries}
      LAPACK_LIBRARIES
    OPTIONAL
      MPI_CXX_LIBRARIES
      MPI_Fortran_LIBRARIES
  INCLUDE_DIRS
    REQUIRED PSBLAS_INCLUDE_DIR
  CLEAR
    ${_psblas_library_variables} ${_additional_libraries} PSBLAS_INCLUDE_DIR
  )
