// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/qr.h>

#include <complex>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace QRImplementation
  {
    // see the corresponding note in the header
    template <typename Number>
    void
    call_trmv(const char            uplo,
              const char            trans,
              const char            diag,
              const types::blas_int n,
              const Number         *a,
              const types::blas_int lda,
              Number               *x,
              const types::blas_int incx)
    {
      trmv(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
    }

    template <typename Number>
    void
    call_trtrs(const char            uplo,
               const char            trans,
               const char            diag,
               const types::blas_int n,
               const types::blas_int nrhs,
               const Number         *a,
               const types::blas_int lda,
               Number               *b,
               const types::blas_int ldb,
               types::blas_int      *info)
    {
      trtrs(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, info);
    }

    template void
    call_trmv(const char,
              const char,
              const char,
              const types::blas_int,
              const float *,
              const types::blas_int,
              float *,
              const types::blas_int);

    template void
    call_trmv(const char,
              const char,
              const char,
              const types::blas_int,
              const double *,
              const types::blas_int,
              double *,
              const types::blas_int);

    template void
    call_trmv(const char,
              const char,
              const char,
              const types::blas_int,
              const std::complex<float> *,
              const types::blas_int,
              std::complex<float> *,
              const types::blas_int);

    template void
    call_trmv(const char,
              const char,
              const char,
              const types::blas_int,
              const std::complex<double> *,
              const types::blas_int,
              std::complex<double> *,
              const types::blas_int);

    template void
    call_trtrs(const char,
               const char,
               const char,
               const types::blas_int,
               const types::blas_int,
               const float *,
               const types::blas_int,
               float *,
               const types::blas_int,
               types::blas_int *);

    template void
    call_trtrs(const char,
               const char,
               const char,
               const types::blas_int,
               const types::blas_int,
               const double *,
               const types::blas_int,
               double *,
               const types::blas_int,
               types::blas_int *);

    template void
    call_trtrs(const char,
               const char,
               const char,
               const types::blas_int,
               const types::blas_int,
               const std::complex<float> *,
               const types::blas_int,
               std::complex<float> *,
               const types::blas_int,
               types::blas_int *);

    template void
    call_trtrs(const char,
               const char,
               const char,
               const types::blas_int,
               const types::blas_int,
               const std::complex<double> *,
               const types::blas_int,
               std::complex<double> *,
               const types::blas_int,
               types::blas_int *);
  } // namespace QRImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
