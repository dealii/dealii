// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
              const Number *        a,
              const types::blas_int lda,
              Number *              x,
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
               const Number *        a,
               const types::blas_int lda,
               Number *              b,
               const types::blas_int ldb,
               types::blas_int *     info)
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
