## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2019 by the deal.II authors
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
# Configuration for the lapack library:
#

MACRO(FEATURE_LAPACK_FIND_EXTERNAL var)
  CLEAR_CMAKE_REQUIRED()
  FIND_PACKAGE(LAPACK)

  #
  # We do a check for availability of every single LAPACK function we use.
  #
  IF(LAPACK_FOUND)
    SET(${var} TRUE)

    #
    # Clear the test flags because the following test will use a C compiler
    #
    CLEAR_CMAKE_REQUIRED()
    SET(CMAKE_REQUIRED_LIBRARIES
      ${DEAL_II_LINKER_FLAGS_SAVED} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
      )

    CHECK_C_SOURCE_COMPILES("
      char daxpy_(); char dgeev_(); char dgeevx_(); char dgelsd_(); char
      dgemm_(); char dgemv_(); char dgeqrf_(); char dgesdd_(); char
      dgesvd_(); char dgetrf_(); char dgetri_(); char dgetrs_(); char
      dorgqr_(); char dormqr_(); char dstev_(); char dsyevx_(); char
      dsygv_(); char dsygvx_(); char dtrtrs_(); char saxpy_(); char
      sgeev_(); char sgeevx_(); char sgelsd_(); char sgemm_(); char
      sgemv_(); char sgeqrf_(); char sgesdd_(); char sgesvd_(); char
      sgetrf_(); char sgetri_(); char sgetrs_(); char sorgqr_(); char
      sormqr_(); char sstev_(); char ssyevx_(); char ssygv_(); char
      ssygvx_(); char strtrs_();
      int main(){
        daxpy_ (); dgeev_ (); dgeevx_ (); dgelsd_ (); dgemm_ (); dgemv_ ();
        dgeqrf_ (); dgesdd_ (); dgesvd_ (); dgetrf_ (); dgetri_ (); dgetrs_
        (); dorgqr_ (); dormqr_ (); dstev_ (); dsyevx_ (); dsygv_ ();
        dsygvx_ (); dtrtrs_ (); saxpy_ (); sgeev_ (); sgeevx_ (); sgelsd_
        (); sgemm_ (); sgemv_ (); sgeqrf_ (); sgesdd_ (); sgesvd_ ();
        sgetrf_ (); sgetri_ (); sgetrs_ (); sorgqr_ (); sormqr_ (); sstev_
        (); ssyevx_ (); ssygv_ (); ssygvx_ (); strtrs_ ();

        return 0;
      }"
      LAPACK_SYMBOL_CHECK)

    IF(NOT LAPACK_SYMBOL_CHECK)
      MESSAGE(STATUS
        "Could not find a sufficient BLAS/LAPACK installation: "
        "BLAS/LAPACK symbol check failed! Consult CMakeFiles/CMakeError.log "
        "for further information."
        )
      SET(LAPACK_ADDITIONAL_ERROR_STRING
        ${LAPACK_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient BLAS/LAPACK installation: \n"
        "BLAS/LAPACK symbol check failed! This usually means that your "
        "BLAS/LAPACK installation is incomplete or the link line is "
        "broken. Consult\n"
        "  CMakeFiles/CMakeError.log\n"
        "for further information.\n"
        )
      SET(${var} FALSE)
    ENDIF()

    #
    # See if we use Intel-MKL by compiling a small program
    #
    SET(CMAKE_REQUIRED_INCLUDES
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
    IF(MKL_SYMBOL_CHECK)
      MESSAGE(STATUS
      "Use Intel MKL for BLAS/LAPACK."
      )
      SET(DEAL_II_LAPACK_WITH_MKL ON)
    ELSE()
      MESSAGE(STATUS
      "Use other than Intel MKL implementation of BLAS/LAPACK (consult CMakeFiles/CMakeError.log for further information)."
      )
      SET(DEAL_II_LAPACK_WITH_MKL OFF)
    ENDIF()

  ENDIF()
ENDMACRO()


CONFIGURE_FEATURE(LAPACK)
