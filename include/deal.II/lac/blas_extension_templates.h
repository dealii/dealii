// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_blas_extension_templates_h
#define dealii_blas_extension_templates_h

#include <deal.II/base/config.h>

#include <deal.II/lac/lapack_support.h>

#ifdef DEAL_II_HAVE_FP_EXCEPTIONS
#  include <cfenv>
#endif

// Intel-MKL specific functions
#ifdef DEAL_II_LAPACK_WITH_MKL
// see
// https://software.intel.com/en-us/mkl-windows-developer-guide-using-complex-types-in-c-c
#  define MKL_Complex8 std::complex<float>
#  define MKL_Complex16 std::complex<double>
#  include <mkl_trans.h>
#endif


DEAL_II_NAMESPACE_OPEN


template <typename number1, typename number2, typename number3>
inline void
omatcopy(char,
         char,
         dealii::types::blas_int,
         dealii::types::blas_int,
         const number1,
         const number2 *,
         dealii::types::blas_int,
         number3 *,
         dealii::types::blas_int)
{
  DEAL_II_NOT_IMPLEMENTED();
}



inline void
omatcopy(char                    ordering,
         char                    trans,
         dealii::types::blas_int rows,
         dealii::types::blas_int cols,
         const float             alpha,
         const float            *A,
         dealii::types::blas_int lda,
         float                  *B,
         dealii::types::blas_int ldb)
{
#ifdef DEAL_II_LAPACK_WITH_MKL
  mkl_somatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
#else
  (void)ordering;
  (void)trans;
  (void)rows;
  (void)cols;
  (void)alpha;
  (void)A;
  (void)lda;
  (void)B;
  (void)ldb;
  Assert(false, LAPACKSupport::ExcMissing("mkl_somatcopy"));
#endif // DEAL_II_LAPACK_WITH_MKL
}



inline void
omatcopy(char                    ordering,
         char                    trans,
         dealii::types::blas_int rows,
         dealii::types::blas_int cols,
         const double            alpha,
         const double           *A,
         dealii::types::blas_int lda,
         double                 *B,
         dealii::types::blas_int ldb)
{
#ifdef DEAL_II_LAPACK_WITH_MKL
  mkl_domatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
#else
  (void)ordering;
  (void)trans;
  (void)rows;
  (void)cols;
  (void)alpha;
  (void)A;
  (void)lda;
  (void)B;
  (void)ldb;
  Assert(false, LAPACKSupport::ExcMissing("mkl_domatcopy"));
#endif // DEAL_II_LAPACK_WITH_MKL
}



inline void
omatcopy(char                       ordering,
         char                       trans,
         dealii::types::blas_int    rows,
         dealii::types::blas_int    cols,
         const std::complex<float>  alpha,
         const std::complex<float> *A,
         dealii::types::blas_int    lda,
         std::complex<float>       *B,
         dealii::types::blas_int    ldb)
{
#ifdef DEAL_II_LAPACK_WITH_MKL
  mkl_comatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
#else
  (void)ordering;
  (void)trans;
  (void)rows;
  (void)cols;
  (void)alpha;
  (void)A;
  (void)lda;
  (void)B;
  (void)ldb;
  Assert(false, LAPACKSupport::ExcMissing("mkl_comatcopy"));
#endif // DEAL_II_LAPACK_WITH_MKL
}



inline void
omatcopy(char                        ordering,
         char                        trans,
         dealii::types::blas_int     rows,
         dealii::types::blas_int     cols,
         const std::complex<double>  alpha,
         const std::complex<double> *A,
         dealii::types::blas_int     lda,
         std::complex<double>       *B,
         dealii::types::blas_int     ldb)
{
#ifdef DEAL_II_LAPACK_WITH_MKL
  mkl_zomatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
#else
  (void)ordering;
  (void)trans;
  (void)rows;
  (void)cols;
  (void)alpha;
  (void)A;
  (void)lda;
  (void)B;
  (void)ldb;
  Assert(false, LAPACKSupport::ExcMissing("mkl_zomatcopy"));
#endif // DEAL_II_LAPACK_WITH_MKL
}

DEAL_II_NAMESPACE_CLOSE

#endif
