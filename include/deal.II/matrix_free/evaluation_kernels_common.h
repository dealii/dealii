// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_evaluation_kernels_common_h
#define dealii_matrix_free_evaluation_kernels_common_h

#include <deal.II/base/config.h>

#include <deal.II/base/utilities.h>

#include <deal.II/matrix_free/fe_evaluation_data.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  template <int dim, int fe_degree, typename Number, bool is_face>
  inline void
  embed_truncated_into_full_tensor_product(
    const unsigned int                      n_components,
    Number                                 *values_dofs,
    const Number                           *values_dofs_actual,
    FEEvaluationData<dim, Number, is_face> &fe_eval)
  {
    const auto &shape_info = fe_eval.get_shape_info();
    const auto &shape_data = shape_info.data.front();

    const std::size_t dofs_per_comp =
      Utilities::pow(shape_data.fe_degree + 1, dim);
    const std::size_t n_dofs_per_comp = shape_info.dofs_per_component_on_cell;
    const int degree = fe_degree != -1 ? fe_degree : shape_data.fe_degree;

    for (unsigned int c = 0; c < n_components; ++c)
      for (int i = 0, count_p = 0, count_q = 0; i < (dim > 2 ? degree + 1 : 1);
           ++i)
        {
          for (int j = 0; j < (dim > 1 ? degree + 1 - i : 1); ++j)
            {
              for (int k = 0; k < degree + 1 - j - i; ++k, ++count_p, ++count_q)
                values_dofs[c * dofs_per_comp + count_q] =
                  values_dofs_actual[c * n_dofs_per_comp + count_p];
              for (int k = degree + 1 - j - i; k < degree + 1; ++k, ++count_q)
                values_dofs[c * dofs_per_comp + count_q] = Number();
            }
          for (int j = degree + 1 - i; j < degree + 1; ++j)
            for (int k = 0; k < degree + 1; ++k, ++count_q)
              values_dofs[c * dofs_per_comp + count_q] = Number();
        }
  }



  template <int dim, int fe_degree, typename Number, bool is_face>
  inline void
  truncate_tensor_product_to_complete_degrees(
    const unsigned int                      n_components,
    Number                                 *values_dofs_actual,
    const Number                           *values_dofs,
    FEEvaluationData<dim, Number, is_face> &fe_eval)
  {
    const auto &shape_info = fe_eval.get_shape_info();
    const auto &shape_data = shape_info.data.front();

    const unsigned int dofs_per_comp =
      Utilities::fixed_power<dim>(shape_data.fe_degree + 1);
    const std::size_t n_dofs_per_comp = shape_info.dofs_per_component_on_cell;
    const int degree = fe_degree != -1 ? fe_degree : shape_data.fe_degree;
    for (unsigned int c = 0; c < n_components; ++c)
      for (int i = 0, count_p = 0, count_q = 0; i < (dim > 2 ? degree + 1 : 1);
           ++i)
        {
          for (int j = 0; j < (dim > 1 ? degree + 1 - i : 1); ++j)
            {
              for (int k = 0; k < degree + 1 - j - i; ++k, ++count_p, ++count_q)
                values_dofs_actual[c * n_dofs_per_comp + count_p] =
                  values_dofs[c * dofs_per_comp + count_q];
              count_q += j + i;
            }
          count_q += i * (degree + 1);
        }
  }
} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
