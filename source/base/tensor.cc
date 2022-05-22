// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2022 by the deal.II authors
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

#include <deal.II/base/tensor.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/lapack_templates.h>

#include <array>

DEAL_II_NAMESPACE_OPEN

namespace
{
  template <int dim, typename Number>
  void
  calculate_svd_in_place(Tensor<2, dim, Number> &A_in_VT_out,
                         Tensor<2, dim, Number> &U)
  {
    // inputs: A
    // outputs: V^T, U
    // SVD: A = U S V^T
    // Since Tensor stores data in row major order and lapack expects column
    // major ordering, we take the SVD of A^T by running the gesvd command.
    // The results (V^T)^T and U^T are provided in column major that we use
    // as row major results V^T and U.
    // It essentially computs A^T = (V^T)^T S U^T and gives us V^T and U.
    // This trick gives what we originally wanted (A = U S V^T) but the order
    // of U and V^T is reversed.
    std::array<Number, dim> S;
    const types::blas_int   N = dim;
    // lwork must be >= max(1, 3*min(m,n)+max(m,n), 5*min(m,n))
    const types::blas_int     lwork = 5 * dim;
    std::array<Number, lwork> work;
    types::blas_int           info;
    gesvd(&LAPACKSupport::O, // replace VT in place
          &LAPACKSupport::A,
          &N,
          &N,
          A_in_VT_out.begin_raw(),
          &N,
          S.data(),
          A_in_VT_out.begin_raw(),
          &N,
          U.begin_raw(),
          &N,
          work.data(),
          &lwork,
          &info);
    Assert(info == 0, LAPACKSupport::ExcErrorCode("gesvd", info));
    Assert(S.back() / S.front() > 1.e-10, LACExceptions::ExcSingular());
  }
} // namespace



template <int dim, typename Number>
Tensor<2, dim, Number>
project_onto_orthogonal_tensors(const Tensor<2, dim, Number> &A)
{
  Tensor<2, dim, Number> VT(A), U;
  calculate_svd_in_place(VT, U);
  return U * VT;
}



template Tensor<2, 1, float>
project_onto_orthogonal_tensors(const Tensor<2, 1, float> &);
template Tensor<2, 2, float>
project_onto_orthogonal_tensors(const Tensor<2, 2, float> &);
template Tensor<2, 3, float>
project_onto_orthogonal_tensors(const Tensor<2, 3, float> &);
template Tensor<2, 1, double>
project_onto_orthogonal_tensors(const Tensor<2, 1, double> &);
template Tensor<2, 2, double>
project_onto_orthogonal_tensors(const Tensor<2, 2, double> &);
template Tensor<2, 3, double>
project_onto_orthogonal_tensors(const Tensor<2, 3, double> &);

DEAL_II_NAMESPACE_CLOSE
