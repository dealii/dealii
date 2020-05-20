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
#include <deal.II/lac/utilities.h>

#include <complex>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace UtilitiesImplementation
  {
    // see the corresponding note in the header
    template <typename Number>
    void
    call_stev(const char            jobz,
              const types::blas_int n,
              Number *              d,
              Number *              e,
              Number *              z,
              const types::blas_int ldz,
              Number *              work,
              types::blas_int *     info)
    {
      stev(&jobz, &n, d, e, z, &ldz, work, info);
    }


    template void
    call_stev(const char,
              const types::blas_int,
              float *,
              float *,
              float *,
              const types::blas_int,
              float *,
              types::blas_int *);

    template void
    call_stev(const char,
              const types::blas_int,
              double *,
              double *,
              double *,
              const types::blas_int,
              double *,
              types::blas_int *);

    template void
    call_stev(const char,
              const types::blas_int,
              std::complex<float> *,
              std::complex<float> *,
              std::complex<float> *,
              const types::blas_int,
              std::complex<float> *,
              types::blas_int *);

    template void
    call_stev(const char,
              const types::blas_int,
              std::complex<double> *,
              std::complex<double> *,
              std::complex<double> *,
              const types::blas_int,
              std::complex<double> *,
              types::blas_int *);
  } // namespace UtilitiesImplementation
} // namespace internal



DEAL_II_NAMESPACE_CLOSE
