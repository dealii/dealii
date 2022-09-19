// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

#ifndef dealii_tensor_product_matrix_templates_h
#define dealii_tensor_product_matrix_templates_h


#include <deal.II/base/config.h>

#include <deal.II/lac/tensor_product_matrix.h>

#ifndef FDM_DEGREE_MAX
// We set this value 17, since this is the value needed for smoothers for
// cell-centered patches with overlap for continuous elements with degrees up
// to FE_EVAL_FACTORY_DEGREE_MAX.
#  define FDM_DEGREE_MAX 17
#endif

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace TensorProductMatrixSymmetricSum
  {
    template <int n_rows_1d_templated, std::size_t dim, typename Number>
    void
    select_vmult(Number *                               dst,
                 const Number *                         src,
                 AlignedVector<Number> &                tmp,
                 const unsigned int                     n_rows_1d,
                 const std::array<const Number *, dim> &mass_matrix,
                 const std::array<const Number *, dim> &derivative_matrix)
    {
      if (n_rows_1d_templated == n_rows_1d)
        vmult<n_rows_1d_templated>(
          dst, src, tmp, n_rows_1d, mass_matrix, derivative_matrix);
      else if (n_rows_1d_templated < FDM_DEGREE_MAX)
        select_vmult<std::min(n_rows_1d_templated + 1, FDM_DEGREE_MAX)>(
          dst, src, tmp, n_rows_1d, mass_matrix, derivative_matrix);
      else
        vmult<0>(dst, src, tmp, n_rows_1d, mass_matrix, derivative_matrix);
    }



    template <int n_rows_1d_templated, std::size_t dim, typename Number>
    void
    select_apply_inverse(Number *                               dst,
                         const Number *                         src,
                         AlignedVector<Number> &                tmp,
                         const unsigned int                     n_rows_1d,
                         const std::array<const Number *, dim> &eigenvectors,
                         const std::array<const Number *, dim> &eigenvalues)
    {
      if (n_rows_1d_templated == n_rows_1d)
        apply_inverse<n_rows_1d_templated>(
          dst, src, tmp, n_rows_1d, eigenvectors, eigenvalues);
      else if (n_rows_1d_templated < FDM_DEGREE_MAX)
        select_apply_inverse<std::min(n_rows_1d_templated + 1, FDM_DEGREE_MAX)>(
          dst, src, tmp, n_rows_1d, eigenvectors, eigenvalues);
      else
        apply_inverse<0>(dst, src, tmp, n_rows_1d, eigenvectors, eigenvalues);
    }
  } // namespace TensorProductMatrixSymmetricSum
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
