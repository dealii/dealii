## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
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
# Configuration for the lapack library:
#

macro(feature_lapack_find_external var)
  clear_cmake_required()
  find_package(DEAL_II_LAPACK)

  #
  # We do a check for availability of every single LAPACK function we use.
  #
  if(LAPACK_FOUND)
    set(${var} TRUE)

    #
    # Clear the test flags
    #
    clear_cmake_required()
    set(CMAKE_REQUIRED_LIBRARIES
      ${DEAL_II_LINKER_FLAGS_SAVED} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
      )

    CHECK_CXX_SOURCE_COMPILES("
      #define MANGLE(name, NAME) ${DEAL_II_FORTRAN_MANGLE}

      // Declare BLAS and LAPACK functions
      extern \"C\" {
        char MANGLE(daxpy, DAXPY)(); char MANGLE(dgeev, DGEEV)(); char MANGLE(dgeevx, DGEEVX)(); char MANGLE(dgelsd, DGELSD)();
        char MANGLE(dgemm, DGEMM)(); char MANGLE(dgemv, DGEMV)(); char MANGLE(dgeqrf, DGEQRF)(); char MANGLE(dgesdd, DGESDD)();
        char MANGLE(dgesvd, DGESVD)(); char MANGLE(dgetrf, DGETRF)(); char MANGLE(dgetri, DGETRI)(); char MANGLE(dgetrs, DGETRS)();
        char MANGLE(dorgqr, DORGQR)(); char MANGLE(dormqr, DORMQR)(); char MANGLE(dstev, DSTEV)(); char MANGLE(dsyevx, DSYEVX)();
        char MANGLE(dsygv, DSYGV)(); char MANGLE(dsygvx, DSYGVX)(); char MANGLE(dtrtrs, DTRTRS)(); char MANGLE(saxpy, SAXPY)();
        char MANGLE(sgeev, SGEEV)(); char MANGLE(sgeevx, SGEEVX)(); char MANGLE(sgelsd, SGELSD)(); char MANGLE(sgemm, SGEMM)();
        char MANGLE(sgemv, SGEMV)(); char MANGLE(sgeqrf, SGEQRF)(); char MANGLE(sgesdd, SGESDD)(); char MANGLE(sgesvd, SGESVD)();
        char MANGLE(sgetrf, SGETRF)(); char MANGLE(sgetri, SGETRI)(); char MANGLE(sgetrs, SGETRS)(); char MANGLE(sorgqr, SORGQR)();
        char MANGLE(sormqr, SORMQR)(); char MANGLE(sstev, SSTEV)(); char MANGLE(ssyevx, SSYEVX)(); char MANGLE(ssygv, SSYGV)();
        char MANGLE(ssygvx, SSYGVX)(); char MANGLE(strtrs, STRTRS)();
      }

      // Then call them. When we link this program, the linker needs to be able
      // to resolve these function references.
      int main(){
        MANGLE(daxpy, DAXPY)(); MANGLE(dgeev, DGEEV)(); MANGLE(dgeevx, DGEEVX)(); MANGLE(dgelsd, DGELSD)(); MANGLE(dgemm, DGEMM)();
        MANGLE(dgemv, DGEMV)(); MANGLE(dgeqrf, DGEQRF)(); MANGLE(dgesdd, DGESDD)(); MANGLE(dgesvd, DGESVD)(); MANGLE(dgetrf, DGETRF)();
        MANGLE(dgetri, DGETRI)(); MANGLE(dgetrs, DGETRS)(); MANGLE(dorgqr, DORGQR)(); MANGLE(dormqr, DORMQR)(); MANGLE(dstev, DSTEV)();
        MANGLE(dsyevx, DSYEVX)(); MANGLE(dsygv, DSYGV)(); MANGLE(dsygvx, DSYGVX)(); MANGLE(dtrtrs, DTRTRS)(); MANGLE(saxpy, SAXPY)();
        MANGLE(sgeev, SGEEV)(); MANGLE(sgeevx, SGEEVX)(); MANGLE(sgelsd, SGELSD)(); MANGLE(sgemm, SGEMM)(); MANGLE(sgemv, SGEMV)();
        MANGLE(sgeqrf, SGEQRF)(); MANGLE(sgesdd, SGESDD)(); MANGLE(sgesvd, SGESVD)(); MANGLE(sgetrf, SGETRF)(); MANGLE(sgetri, SGETRI)();
        MANGLE(sgetrs, SGETRS)(); MANGLE(sorgqr, SORGQR)(); MANGLE(sormqr, SORMQR)(); MANGLE(sstev, SSTEV)(); MANGLE(ssyevx, SSYEVX)();
        MANGLE(ssygv, SSYGV)(); MANGLE(ssygvx, SSYGVX)(); MANGLE(strtrs, STRTRS)();

        return 0;
      }"
      LAPACK_SYMBOL_CHECK)

    if(NOT LAPACK_SYMBOL_CHECK)
      message(STATUS
        "Could not find a sufficient BLAS/LAPACK installation: "
        "BLAS/LAPACK symbol check failed! Consult CMakeFiles/CMakeError.log "
        "for further information."
        )
      set(LAPACK_ADDITIONAL_ERROR_STRING
        ${LAPACK_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient BLAS/LAPACK installation: \n"
        "BLAS/LAPACK symbol check failed! This usually means that your "
        "BLAS/LAPACK installation is incomplete or the link line is "
        "broken. Consult\n"
        "  CMakeFiles/CMakeError.log\n"
        "for further information.\n"
        )
      set(${var} FALSE)
    endif()

    #
    # See if we use Intel-MKL by compiling a small program
    #
    set(CMAKE_REQUIRED_INCLUDES
      ${LAPACK_INCLUDE_DIRS}
      )
    CHECK_CXX_SOURCE_COMPILES("
    #include <mkl.h>
    #include <vector>
    int main(){
      const int m = 5;
      const int n = 2;
      std::vector<double> A(m*n,0.);
      std::vector<double> B(m*n,0.);
      mkl_domatcopy('C', 'T', m, n, 1., A.data(), n, B.data(), m);
      return 0;
    }"
    MKL_SYMBOL_CHECK)
    if(MKL_SYMBOL_CHECK)
      message(STATUS
      "Use Intel MKL for BLAS/LAPACK."
      )
      set(DEAL_II_LAPACK_WITH_MKL ON)
    else()
      message(STATUS
      "Use other than Intel MKL implementation of BLAS/LAPACK (consult CMakeFiles/CMakeError.log for further information)."
      )
      set(DEAL_II_LAPACK_WITH_MKL OFF)
    endif()

  endif()
endmacro()


configure_feature(LAPACK)
