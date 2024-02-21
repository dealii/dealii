// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tensor_product_matrix_templates_h
#define dealii_tensor_product_matrix_templates_h


#include <deal.II/base/config.h>

#include <deal.II/lac/tensor_product_matrix.h>

#ifndef FDM_N_ROWS_MAX
// Maximum number of rows with pre-compiled TensorProductMatrixSymmetricSum
// kernels. If no value is given by the user during
// compilation, we choose its value so that all number of rows are pre-compiled
// to support smoothers for cell-centered patches with overlap for continuous
// elements with degrees up to FE_EVAL_FACTORY_DEGREE_MAX (default value 6).
#  ifndef FE_EVAL_FACTORY_DEGREE_MAX
#    define FDM_N_ROWS_MAX 17
#  else
#    define FDM_N_ROWS_MAX (FE_EVAL_FACTORY_DEGREE_MAX * 3 - 1)
#  endif
#endif

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace TensorProductMatrixSymmetricSum
  {
    template <int n_rows_1d_templated, std::size_t dim, typename Number>
    void
    select_vmult(Number                                *dst,
                 const Number                          *src,
                 AlignedVector<Number>                 &tmp,
                 const unsigned int                     n_rows_1d,
                 const std::array<const Number *, dim> &mass_matrix,
                 const std::array<const Number *, dim> &derivative_matrix)
    {
      if (n_rows_1d_templated == n_rows_1d)
        vmult<n_rows_1d_templated>(
          dst, src, tmp, n_rows_1d, mass_matrix, derivative_matrix);
      else if (n_rows_1d_templated < FDM_N_ROWS_MAX)
        select_vmult<std::min(n_rows_1d_templated + 1, FDM_N_ROWS_MAX)>(
          dst, src, tmp, n_rows_1d, mass_matrix, derivative_matrix);
      else
        vmult<0>(dst, src, tmp, n_rows_1d, mass_matrix, derivative_matrix);
    }



    template <int n_rows_1d_templated, std::size_t dim, typename Number>
    void
    select_apply_inverse(Number                                *dst,
                         const Number                          *src,
                         const unsigned int                     n_rows_1d,
                         const std::array<const Number *, dim> &eigenvectors,
                         const std::array<const Number *, dim> &eigenvalues,
                         const Number *inverted_eigenvalues)
    {
      if (n_rows_1d_templated == n_rows_1d)
        apply_inverse<n_rows_1d_templated>(
          dst, src, n_rows_1d, eigenvectors, eigenvalues, inverted_eigenvalues);
      else if (n_rows_1d_templated < FDM_N_ROWS_MAX)
        select_apply_inverse<std::min(n_rows_1d_templated + 1, FDM_N_ROWS_MAX)>(
          dst, src, n_rows_1d, eigenvectors, eigenvalues, inverted_eigenvalues);
      else
        apply_inverse<0>(
          dst, src, n_rows_1d, eigenvectors, eigenvalues, inverted_eigenvalues);
    }
  } // namespace TensorProductMatrixSymmetricSum
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
